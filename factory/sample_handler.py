import pickle
import os
import re

import hist
import numpy as np
import scipy.integrate as integrate
from scipy.stats.distributions import chi2
import uproot
import json

from config import Configurator
from ntuple_loader import NTupleLoader
import utils
from utils import Ellipse, REGION_ID_MAP, ELLIPSE_X, ELLIPSE_Y, CMSSW_BASE, TEMPLATES


class SampleSet():
    """
    A SampleSet corresponds to one output file. There
    is one SampleSet for each unique combination of
    [Signal, Channel, Region, Sensitive Variable]
    """

    def __init__(self, set_params: dict):
        print(f"Initializing {set_params['signal']} "
              f"{set_params['channel']} {set_params['region']} "
              f"{set_params['year']} {set_params['svar']}")
        self.set_params = set_params
        self.config = Configurator()
        self.samples = {
            "signal": Sample(set_params["signal"], set_params),
            "data": Sample("data", set_params),
        } | {bkg: Sample(bkg, set_params) for bkg in self.config.backgrounds}
        self.set_binning = None

    def create_hists_and_save_file(self):
        self._create_histograms()
        self._assert_nonzero_backgrounds()
        self._save_output_file()

    def set_sample_bins(self, ellipses: Ellipse = None):
        """
        Sets the binning to the binning obtained from
        the binner in all samples in the set.
        """
        if ellipses:
            self.set_binning = ellipses
        else:
            binner = Binner(self.samples["signal"], self.set_params)
            self.set_binning = binner.get_binning()

        for sample in self.samples.values():
            sample.bins = self.set_binning

    def _create_histograms(self):
        """
        Calls all the samples to create all relevant
        histograms for the output root file.
        """
        for sample in self.samples.values():
            sample.create_histograms()

    def _merge_bins_of_sample(self, hist_tuple, idx_zero_bin, idx_neighbor_bin, idx_merged_bin):
        """
        Arguments:
        hist_tuple: A tuple of (bin_yield, bin_yield_sq_weights,
             bin_edges)
        idx_zero_bin: The index of the bin with 0 yield that is
            to be merged with a neighboring bin
        idx_neighbor_bin: The index of the bin that is to be
            merged with the bin at idx_zero_bin
        idx_merged_bin: The index of the merged bin in the
            array resulting from the merge.

        Returns:
        A modified version of the input hist_tuple with the one less
        bin, due to the bin at index idx_zero_bin being merged with
        its neighbor.
        """
        old_bin_content = hist_tuple[0]
        old_bin_content_sq_weights = hist_tuple[1]
        old_bin_boundries = hist_tuple[2]
        new_bin_values = np.delete(old_bin_content, idx_zero_bin)
        new_bin_values_sq_weights = np.delete(old_bin_content_sq_weights, idx_zero_bin)
        if old_bin_content[idx_zero_bin] < 0:
            new_bin_values[idx_merged_bin] += old_bin_content[idx_zero_bin]
            new_bin_values_sq_weights[idx_merged_bin] += old_bin_content_sq_weights[idx_zero_bin]
        new_bin_boundries = np.delete(old_bin_boundries, idx_neighbor_bin)
        return (new_bin_values, new_bin_values_sq_weights, new_bin_boundries)

    def _merge_bins(self, bkg_vals):
        """
        This function recursively merges the bins that
        have 0 content in the sum of the nominal background
        samples.
        """
        # If there are no 0 bkg bins; do nothing!
        if (bkg_vals > 0).all():
            return bkg_vals

        # Find indeces of bins before and after merge
        n_bins = bkg_vals.shape[0]
        idx_zero_bin = np.where(bkg_vals <= 0)[0][0]
        idx_neighbor_bin = idx_zero_bin + 1 if n_bins > idx_zero_bin + 1 else idx_zero_bin - 1
        idx_merged_bin = idx_zero_bin if n_bins > idx_zero_bin + 1 else idx_zero_bin - 1

        new_nominal_values = np.zeros_like(bkg_vals[1:])
        for sample_name, sample in self.samples.items():
            # Merge bins of variation histograms
            for var_name, var in sample.var_hists.items():
                self.samples[sample_name].var_hists[var_name] = self._merge_bins_of_sample(
                    var, idx_zero_bin, idx_neighbor_bin, idx_merged_bin
                )

            # Merge bins of nominal histograms
            self.samples[sample_name].svar_hist = self._merge_bins_of_sample(
                sample.svar_hist, idx_zero_bin, idx_neighbor_bin, idx_merged_bin
            )
            if sample_name not in ["data", "signal"]:
                new_nominal_values += self.samples[sample_name].svar_hist[0]

        if (new_nominal_values >= 0).all():
            assert np.sum(new_nominal_values) == np.sum(bkg_vals)

        if (new_nominal_values > 0).all():
            print("merged_two_bin")
            print(self.set_params)
            return new_nominal_values

        return self._merge_bins(new_nominal_values)

    def _assert_nonzero_backgrounds(self):
        """
        A validation function that checks after all the hists
        in the sample set are built, that their sum does not
        have empty bins in the signal region.
        """
        bkg_vals = np.zeros_like(self.samples["data"].svar_hist[0])
        for sname, bkg_s in self.samples.items():
            if sname in ["data", "signal"]:
                continue
            bkg_vals += bkg_s.svar_hist[0]

        is_signal_region = "SR" in self.set_params["region"]
        if not (bkg_vals > 0).all() and not is_signal_region:
            bkg_vals = self._merge_bins(bkg_vals)

        assert (bkg_vals > 0).all(), (
            f"Bin with no BKG for {self.set_params['svar']} in "
            f"{self.set_params['region']} region and {self.set_params['channel']}."
        )

    def _ntuple_var_to_combine_key(self, v: str, sname: str):
        """
        The keys that UHH2 uses to store the weights
        for the variations are not the same we want to
        use for the variation histograms in the datacards,
        they can however be maped to the combine keys,
        which is what this function does.
        """
        # Remove "weight" prefixes
        v = re.sub(r"^weight_", "", v)
        # Format postfixes correctly
        v = re.sub(r"pdf_up$", f"pdf_{sname}Up", v)
        v = re.sub(r"pdf_down$", f"pdf_{sname}Down", v)
        v = re.sub(r"pu_(\w+)", r"pileup_YEAR_\1", v)
        v = re.sub(r"puid_sf_(\w+)", r"pileupJetID_YEAR_\1", v)
        v = re.sub(r"btag_bc_(\w+)$", r"btag_bc_YEAR_\1", v)
        v = re.sub(r"btag_light_(\w+)$", r"btag_light_YEAR_\1", v)
        v = re.sub(r"_up$", "Up", v)
        v = re.sub(r"_d(ow)?n$", "Down", v)
        v = re.sub(r"muf_noneup$", f"muf_{sname}Up", v)
        v = re.sub(r"muf_nonedown$", f"muf_{sname}Down", v)
        v = re.sub(r"mur_upnone$", f"mur_{sname}Up", v)
        v = re.sub(r"mur_downnone$", f"mur_{sname}Down", v)
        v = re.sub(r"murmuf_upup$", f"murmuf_{sname}Up", v)
        v = re.sub(r"murmuf_downdown$", f"murmuf_{sname}Down", v)
        # Insert "CMS_"
        if "fsr" not in v and "isr" not in v and "pdf" not in v:
            v = re.sub(r"^([a-z]+_)", r"CMS_\1", v)
        v = re.sub(r"^(pileup)", r"CMS_\1", v)
        v = re.sub(r"^(btag)", r"CMS_\1", v)
        v = re.sub(r"^(btag_bc)", r"CMS_\1", v)
        v = re.sub(r"^(btag_light)", r"CMS_\1", v)
        v = re.sub(r"^(btag)_central", r"CMS_\1", v)
        v = re.sub(r"^(jer)", r"CMS_res_j_YEAR", v)
        v = re.sub(r"^(mur)$", r"CMS_\1", v)
        v = re.sub(r"^(muf)$", r"CMS_\1", v)
        v = re.sub(r"^(murmuf)", r"CMS_\1", v)
        # JEC
        v = re.sub(r"(YEAR)", rf"{self.set_params['year']}", v)
        v = re.sub(r"^(Total)", r"CMS_scale_j_Total", v)
        v = re.sub(r"^(Absolute)", r"CMS_scale_j_Absolute", v)
        v = re.sub(r"^(BBEC1)", r"CMS_scale_j_BBEC1", v)
        v = re.sub(r"^(EC2)", r"CMS_scale_j_EC2", v)
        v = re.sub(r"^(FlavorQCD)", r"CMS_scale_j_FlavorQCD", v)
        v = re.sub(r"^(HF)", r"CMS_scale_j_HF", v)
        v = re.sub(r"^(RelativeBal)", r"CMS_scale_j_RelativeBal", v)
        v = re.sub(r"^(RelativeSample)", r"CMS_scale_j_RelativeSample", v)

        return v

    def _build_hist_object(self, hist_tuple):
        """
        Builds a hist.Hist type histogram with the
        correct statistical errors. This is necessary
        since Combine uses the TH1 method `GetBinError`
        to get the statistical error for the fits. Simply
        writing a NumPy histogram via uproot will yield
        errant results.

        Arguments:
        hist_tuple[0] -- bin content of the histogram (yield)
        hist_tuple[1] -- bin content of the histogram built with
            weights=weights**2. This is equivalent to
            the per bin variance.
        hist_tuple[2] -- bin boundries of the histogram

        Returns:
        hist.Hist type histogram with correct
        """
        h_yield, h_yield_sq_weights, bin_edges = hist_tuple
        #if "overflow" in bin_edges: bin_edges = [i for i in range(len(bin_edges))]
        h = hist.Hist.new.Var(bin_edges, name="", label="")
        h = h.Weight()
        h.values()[:] = h_yield
        h.variances()[:] = h_yield_sq_weights
        return h

    def _write_var_hists_to_file(self, f, var_hists, region, sname):
        variations = [
            x + xvar for x, xvars in self.config.variations.items()
            for xvar in xvars if xvar and "central" not in xvar and "vjets" not in x
        ]

        vjet_vars = [
            x + xvar for x, xvars in self.config.variations.items()
            for xvar in xvars if "vjets" in x and xvar != ""
        ]
        vjet_vars = [var + _dir for var in vjet_vars for _dir in ("_up", "_dn")]
        variations = variations + vjet_vars

        if self.config.sample_variations:
            relevant_vars = {
                x: y for x, y in
                self.config.sample_variations.items() if
                sname in self.config.sample_var_whitelist[x]
            }
            sample_variations = [
                x + xvar for x, xvars in relevant_vars.items()
                for xvar in xvars if xvar
            ]
        else:
            sample_variations = []

        for v in variations + sample_variations:
            combine_key = self._ntuple_var_to_combine_key(v, sname)
            f[f"{region}/{sname}_{combine_key}"] = self._build_hist_object(var_hists[v])
            if v in sample_variations:
                f[f"{region}/{sname}_{combine_key}_unweighted"] = self._build_hist_object(var_hists[f"{v}_unweighted"])

    def _save_output_file(self):
        """
        Writes the histograms from the samples
        in the set to a .root file with the
        appropriate structure.
        """
        print(f"Saving output ... {self.set_params['signal']} "
              f"{self.set_params['channel']} {self.set_params['region']} "
              f"{self.set_params['year']} {self.set_params['svar']}")
        sgnl = self.set_params["signal"]
        channel = self.set_params["channel"]
        region = self.set_params["region"]
        svar = utils.OUTPUT_MAP[self.set_params["svar"]]

        base_path = os.path.join(TEMPLATES, self.set_params["year"])
        fpath_met = os.path.join(base_path, f"{sgnl}_{svar}_{channel}_{region}.root")

        with uproot.recreate(fpath_met) as f:
            for sname, sample in self.samples.items():
                if sname in ["data", "signal"]:
                    continue
                if np.sum(sample.svar_hist[0]) > 0.5:
                    f[region + '/' + sname] = self._build_hist_object(sample.svar_hist)
                    self._write_var_hists_to_file(f, sample.var_hists, region, sname)
            f[region + "/data_obs"] = self._build_hist_object(self.samples["data"].svar_hist)
            f[region + "/AtoZH"] = self._build_hist_object(self.samples["signal"].svar_hist)
            self._write_var_hists_to_file(f, self.samples["signal"].var_hists, region, "AtoZH")


class Sample():
    """
    A Sample corresponds to a MC physics prcess (e.g. TT)
    or data. All samples are part of a sample set that
    dicatates which sensitive variable, channel, region
    and year are considered.
    """

    def __init__(self, sample, sample_params):
        self.config = Configurator()
        self.sample_params = sample_params
        sample_loader = NTupleLoader()

        # Sample Parameter
        self.sample = sample
        self.svar = sample_params["svar"]
        self.signal = sample_params["signal"]
        self.channel = sample_params["channel"]
        self.region = sample_params["region"]
        self.year = sample_params["year"]

        # Quick Access to Sample Loader
        self.tree = sample_loader.nominal_trees[self.year][sample]
        self.tree_signal = sample_loader.nominal_trees[self.year][self.signal]
        self.weights = self.tree["event_weight"]

        # Histogramed Sample
        self.bins = None
        self.sel = None
        self.svar_hist = None
        self.var_sel = {}
        self.var_hists = {}

        self._set_svar_selection_mask()

        if sample != "data" and sample_loader.sample_vars[self.year]:
            self.sample_vars = sample_loader.sample_vars[self.year][sample]
            self.sample_vars_signal = sample_loader.sample_vars[self.year][self.signal]
            self._set_sample_vars_selection_masks()

    def _cut_on_angle(self, signal_tree, sample_tree):
        """
        Cut out upper part of [dR(A,H), max(dR(Z,t), dR(Z,tbar))] plane.
        Use a linear function of type `y = slope * x + offset` to split
        the plane. Scale down offset to retain 98% efficiency on signal.
        Returns the selection mask.
        """
        def signal_efficiency(y_offset):
            # TODO: Compute Efficiency only on SR
            sel = slope * signal_tree["angle_ah"] + y_offset > signal_tree["angle_max_tz"]
            return np.sum(sel) / signal_tree["angle_ah"].shape[0]

        slope = 0.8125
        target_signal_efficiency = self.config.angle_cut_signal_efficiency
        y_offset = -4  # initialization, optimized below

        # Sacle offset
        while signal_efficiency(y_offset) < target_signal_efficiency:
            y_offset += .04

        selection_mask = slope * sample_tree["angle_ah"] + y_offset > sample_tree["angle_max_tz"]

        return selection_mask

    def _construct_feature_vector(self, sample_tree, mA, mH):
        a_mass = sample_tree["a_mass_chi_sq_fixedb"]
        h_mass = sample_tree["h_mass_reco_chi_sq_fixedb"]
        z_mass_reco = sample_tree["z_mass_reco"]
        delta_m = sample_tree["a_minus_h_mass_chi_sq_fixedb"]
        chi_sq = sample_tree["chi_sq_lowest"]
        a_mass_param = np.ones_like(a_mass) * mA
        h_mass_param = np.ones_like(a_mass) * mH
        pt_edge = utils.pt_of_z_edge(a_mass_param, h_mass_param)

        all_arrays = [a_mass, h_mass, z_mass_reco, delta_m, chi_sq, a_mass_param, h_mass_param, pt_edge]

        return np.dstack(all_arrays)

    def _cut_on_dnn(self, signal_tree, sample_tree):
        # construct feature vector
        mA = int(re.findall(r"m[AH](\d+)", self.signal)[0])
        mH = int(re.findall(r"m[AH](\d+)", self.signal)[1])
        features = self._construct_feature_vector(sample_tree, mA, mH)
        leading_jet_phis = sample_tree["jetsAk4CHS/jetsAk4CHS.m_phi_1"]

        # call inference tool
        inference_tool = InferenceTool()
        sel = inference_tool.predict(features, leading_jet_phis, mA, mH)

        return sel

    def _set_svar_selection_mask(self):
        """
        self.tree["region"] -> region of all events in analysistree
        self.sample_params["region"] -> region currently computed
        """
        is_region = self.tree["region"] == REGION_ID_MAP[self.sample_params["region"]]
        self.sel = is_region

    def _set_sample_vars_selection_masks(self):
        """
        Sets self.var_sel, which is the mask that is used
        to filter the jet var trees for channel, region and
        if specified fullfilment of angle requirements.
        """
        relevant_vars = {
            x: y for x, y in
            self.config.sample_variations.items() if
            self.sample in self.config.sample_var_whitelist[x]
        }
        sample_variations = [x + xvar for x, xvars in relevant_vars.items() for xvar in xvars]
        for variation in sample_variations:
            # Region Selection
            is_region = self.sample_vars[variation]["region"] == REGION_ID_MAP[self.sample_params["region"]]
            self.var_sel[variation] = is_region


    def _get_raw_twod_inputs(self):
        x = self.tree[ELLIPSE_X][self.sel]
        y = self.tree[ELLIPSE_Y][self.sel]
        return [x, y]

    def _build_histogram(self, inputs, weights, check_weights=False):
        if self.svar == "ellipses":
            return self._build_elliptically_binned_hist_fast(inputs[0], inputs[1], weights)
        else:
            H, bins = np.histogram(inputs, bins=self.bins, weights=weights)
            H_sq, _ = np.histogram(inputs, bins=self.bins, weights=weights * weights)
            return (H, H_sq, bins)

    def _build_svar_hist(self):
        svar_weights = self.weights[self.sel]

        if self.svar == "ellipses":
            self.svar_hist = self._build_histogram(self._get_raw_twod_inputs(), svar_weights, check_weights=True)
        else:
            self.svar_hist = self._build_histogram(self.tree[self.svar][self.sel], weights=svar_weights)

    def _build_rms_hists(self):
        nominal, nominal_sq_weights, bins = self.svar_hist

        if self.svar == "ellipses":
            inputs = self._get_raw_twod_inputs()
        else:
            inputs = self.tree[self.svar][self.sel]

        with open(f"norm_facts/{self.sample}_{self.year}_norm.json") as f:
            data = f.read()
        norm_facts = json.loads(data)
        norm_fact_pdf = norm_facts["PDF"]

        sum_bins = 0.0
        for i in range(100):
            pdf_wi = np.array([item[i] for item in self.tree["weight_pdf"][self.sel]])
            var, _, _ = self._build_histogram(inputs, weights=self.weights[self.sel] * pdf_wi)
            bin_content = var * norm_fact_pdf[i]
            sum_bins += (bin_content - nominal)**2
        rms = np.sqrt(sum_bins / 100)
        self.var_hists["weight_pdf_up"] = (np.array(nominal + rms), np.array(nominal + rms), bins)
        self.var_hists["weight_pdf_down"] = (np.array(nominal - rms), np.array(nominal - rms), bins)

    def _build_vjet_var_hists(self, inputs, nom_key, xvars):
        for xvar in xvars:
            w_nom = np.array(self.tree[nom_key][self.sel])
            dkappa = np.array(self.tree[nom_key + xvar][self.sel])
            w_up = self.weights[self.sel] / w_nom * (w_nom + dkappa)
            w_dn = self.weights[self.sel] / w_nom * (w_nom - dkappa)
            self.var_hists[nom_key + xvar + "_up"] = self._build_histogram(inputs, weights=w_up)
            self.var_hists[nom_key + xvar + "_dn"] = self._build_histogram(inputs, weights=w_dn)

    def _build_var_hists(self):
        if self.svar == "ellipses":
            inputs = self._get_raw_twod_inputs()
        else:
            inputs = self.tree[self.svar][self.sel]

        for x, xvars in self.config.variations.items():
            if 'pdf' in x:
                continue  # handled in self._build_rms_hists()
            nominal_correction = not xvars[0] or xvars[0].endswith("central")
            if nominal_correction:
                nom_key = x if "btag" not in x else "weight_btag_central"
                w_nom = self.tree[nom_key][self.sel]
                xvars = xvars[1:]
            else:
                w_nom = 1

            if "vjets" in x:
                self._build_vjet_var_hists(inputs, nom_key, xvars)
                continue
            else:
                w_up = self.weights[self.sel] / w_nom * self.tree[x + xvars[0]][self.sel]
                w_dn = self.weights[self.sel] / w_nom * self.tree[x + xvars[1]][self.sel]

            if 'mur' in x or 'muf' in x:
                with open(f"norm_facts/{self.sample}_{self.year}_norm.json") as f:
                    data = f.read()
                norm_facts = json.loads(data)
                norm_w_up = norm_facts["QCD"][x[7:] + xvars[0]]
                norm_w_dn = norm_facts["QCD"][x[7:] + xvars[1]]
                w_up *= norm_w_up
                w_dn *= norm_w_dn

            self.var_hists[x + xvars[0]] = self._build_histogram(inputs, weights=w_up)
            self.var_hists[x + xvars[1]] = self._build_histogram(inputs, weights=w_dn)

    def _build_sample_var_hists(self):

        if self.config.sample_variations:
            relevant_vars = {
                x: y for x, y in
                self.config.sample_variations.items() if
                self.sample in self.config.sample_var_whitelist[x]
            }
            sample_variations = [x + xvar for x, xvars in relevant_vars.items() for xvar in xvars]
        else:
            sample_variations = []

        for variation in sample_variations:
            sel = self.var_sel[variation]

            if self.svar == "ellipses":
                x = self.sample_vars[variation][ELLIPSE_X][sel]
                y = self.sample_vars[variation][ELLIPSE_Y][sel]
                inputs =  [x, y]
            else:
                inputs = self.sample_vars[variation][self.svar][sel]

            w_var = self.sample_vars[variation]["event_weight"][sel]
            w_ones = np.ones_like(w_var)
            self.var_hists[variation] = self._build_histogram(inputs, weights=w_var)
            self.var_hists[variation + "_unweighted"] = self._build_histogram(inputs, weights=w_ones)

    def _assert_percentage_of_points_in_bin(self, bin_values):
        weight_sum = np.sum(bin_values)

        pct_in_bin = 0
        stddevs = [x.n_std for x in self.bins]
        for n_std, bin_val in zip(stddevs, bin_values[:-1]):
            pct_in_bin += bin_val / weight_sum
            correct_pct, _ = integrate.quad(utils.normal_distribution, -n_std, n_std)
            assert (x := abs(correct_pct - pct_in_bin)) < .01, \
                f"Bin content off by {round(x*100, 2)}% from {n_std} sigma"

    def _build_elliptically_binned_hist_fast(self, x, y, weights=None):
        """
        Builds a np.histogram with bins corresponding to the
        specified ellipses.
        """
        assert isinstance(weights, np.ndarray)
        bin_values = [0 for _ in self.bins]
        bin_values_sq = [0 for _ in self.bins]
        for i, ell in enumerate(self.bins):
            sel_in = ell.is_point_included(x, y)
            bin_values[i] = np.sum(weights[sel_in])
            bin_values_sq[i] = np.sum(weights[sel_in] ** 2)
            x = x[~sel_in]
            y = y[~sel_in]
            weights = weights[~sel_in]

        return (
            np.array(bin_values),
            np.array(bin_values_sq),
            np.array([0] + [e.n_std for e in self.bins])
        )

    def create_histograms(self):
        self._build_svar_hist()
        if not self.sample == "data":
            if "weight_pdf" in self.config.variations:
                self._build_rms_hists()
            self._build_var_hists()
            self._build_sample_var_hists()


class Binner():
    """
    Implements the logic to find the binnning
    for a SampleSet for a specific comnbination
    of Signal, Sensitive Variable and Region.
    """

    binning_map = {
        "stddevs_ellipses": np.array([0.5, 1, 1.5, 2, 2.5, 3]),
        #"stddevs_ellipses": np.linspace(0.1,3,15),
        "defaults": {
            "z_pt_reco": np.linspace(0, 800, 14),
            "jetsAk4CHS.m_phi_1": np.linspace(-3, 3, 30),
            "jetsAk4CHS.m_pt_1": list(range(50,650+1,50)),
            "jetsAk4CHS.m_pt_2": list(range(50,650+1,50)),
            "jetsAk4CHS.m_pt_3": list(range(50,650+1,50)),
            "jetsAk4CHS.m_pt_4": list(range(50,650+1,50)),
            "jetsAk4CHS.m_eta_1": np.linspace(-2.4, 2.4, 9),
            "slimmedElectronsUSER.m_pt_1": np.linspace(20, 550, 16),
            "slimmedElectronsUSER.m_pt_2": np.linspace(20, 550, 16),
            "slimmedMuonsUSER.m_pt_1": np.linspace(20, 550, 16),
            "slimmedMuonsUSER.m_pt_2": np.linspace(20, 550, 16),
            "z_mass_reco": np.linspace(40, 140, 16),
            "mt_A": np.append(list(range(450,2000+1,50)),10000),
            "mt_H": np.append(list(range(200,1500+1,50)),10000),
            "m_H": np.append(list(range(200,1500+1,50)),10000),
            "MET": np.append(list(range(220,1020+1,40)),10000),
            "HT": np.append(list(range(200,1500+1,50)),10000),
            "score": np.linspace(0,1,10)**0.5,
        },
        "SR_5J": {
            "mt_H": [-2,0]
        },
        "CR": {
            "foo": np.linspace(0, 5000, 2),
            "bar": np.linspace(0, 5000, 2),
        },
    }

    def __init__(self, signal_sample: Sample, set_params: dict):
        self.signal_sample = signal_sample
        self.signal = set_params["signal"]
        self.channel = set_params["channel"]
        self.region = set_params["region"]
        self.year = set_params["year"]
        self.svar = set_params["svar"]
        self.tree = signal_sample.tree
        self.sel = signal_sample.sel

    def _get_binning_oned(self):
        try:
            return self.binning_map[self.region][self.svar]
        except KeyError:
            pass  # Fallback to default for variable
        try:
            return self.binning_map["defaults"][self.svar]
        except KeyError:
            return np.linspace(-10000, 10000, 2)

    def _load_CR_ellipses(self):
        ellipses = []
        for n_std in [1, 1.5, 2, 2.5, 3, "overflow"]:
            with open(f"bkg_ellipses/bkg_ellipse_{n_std}.pickle", "rb") as f:
                ell = pickle.load(f)
            ellipses.append(ell)
        return ellipses

    def _compute_ellipses(self):

        loader = NTupleLoader()
        sel = self.signal_sample.sel
        x = self.tree[ELLIPSE_X][self.sel]
        y = self.tree[ELLIPSE_Y][self.sel]
        weights = loader.nominal_trees[self.year][self.signal]["event_weight"][sel]

        ellipses = []
        for n_std in self.binning_map["stddevs_ellipses"]:
            ell = self._fit_ellipse(x, y, weights, n_std)
            ellipses.append(ell)

        return ellipses

    def get_binning(self):
        if self.svar == "ellipses":
            if "SR" in self.region:
                return self._compute_ellipses()
            else:
                #return self._load_CR_ellipses() #pre-generated ellipses required
                return self._compute_ellipses()
        else:
            return self._get_binning_oned()

    def _fit_ellipse(self, x, y, w, n_std: float):
        w_pos = w
        if min(w) < 0:
            w_pos -= min(w)
        assert x.size == y.size == w.size, "x and y must be the same size"

        # Compute Mean
        mean_x = utils.weighted_median(x, w)
        mean_y = utils.weighted_median(y, w)

        cov = np.cov(x, y, aweights=w_pos)
        vals, vecs = np.linalg.eigh(cov)
        angle = np.degrees(np.arctan2(*vecs[::-1, 0]))
        z = integrate.quad(utils.normal_distribution, -n_std, n_std)[0]
        width, height = 2 * np.sqrt(vals * chi2.ppf(z, df=2))

        ell = Ellipse(mean_x, mean_y, width, height, angle, n_std)
        ell.rescale_to_nstd(x, y, w) 

        return ell