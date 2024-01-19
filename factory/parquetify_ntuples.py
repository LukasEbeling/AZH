#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse
import itertools
import os
import re

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import uproot
import utils

from utils import CACHE, CMSSW_BASE
from config import Configurator




class Parquetifier():

    UHH_OUTPUT_PATH = os.path.join(CMSSW_BASE, "src/UHH2/AZH/data/output_02_reconstruction/")

    def __init__(self, sample=None):
        self.config = Configurator()
        self.ul_years = utils.split_ul16(self.config.years)
        self.samples = [re.sub(r"^MA-", "AZH_", sample)]
        if not sample:
            self.samples = self.config.samples

        self.samples = [s.replace("_ljet","").replace("_bjet","") for s in self.samples]
        self.samples = np.unique(self.samples) #drop duplicates

        self.variables_to_load = self.config.svars + self.config.branches

        os.makedirs(CACHE, exist_ok=True)

        if sample.upper() == "DATA":
            self.load_data()
        else:
            for s in self.samples:
                self.load_trees(s)
                self.load_sample_vars(s)

        if "UL16" in self.config.years:
            self.merge_16_pre_post()

    def _tree_to_np_array(self, f, x):
        # This is the case in which the observable
        # is part of a collection, e.g. jet or lepton pt
        if isinstance(x, list):
            branch = f[f"AnalysisTree/{x[0]}"].array(library="np")
            key = f"{x[0]}_{x[1]}"
            try:
                i_th_of_collection = []
                for y in branch:
                    if y.shape[0] > x[1] - 1:
                        i_th_of_collection.append(y[x[1] - 1])
                    else:
                        i_th_of_collection.append(-1)
                branch = np.array(i_th_of_collection)
            except IndexError:
                pass  # For empty Electorn/Muon collections
        else:
            branch = f[f"AnalysisTree/{x}"].array(library="np")
            key = x
        return key, branch

    def load_data(self):

        for year in self.ul_years:
            print(f"Loading {year} DATA")
            data_path = os.path.join(self.UHH_OUTPUT_PATH, "DATA", year, f"DATA.{year}.root")
            with uproot.open(data_path) as f:
                for x in self.variables_to_load:
                    if "Gen" in x: continue
                    key, branch = self._tree_to_np_array(f, x)
                    pa_table = pa.table({"foo": branch})
                    pq.write_table(pa_table, f"{CACHE}/data_{year}_{key}.parquet")

    def _get_relevant_sample_variations(self, sample):
        relevant_vars = {
            x: y for x, y in
            self.config.sample_variations.items() if
            sample in self.config.sample_var_whitelist[x]
        }
        return [x + xvar for x, xvars in relevant_vars.items() for xvar in xvars]

    def load_sample_vars(self, sample):

        for year in self.ul_years:
            print(f"Loading Jet Variations {year} {sample}")
            sample_variations = self._get_relevant_sample_variations(sample)

            base_path = os.path.join(self.UHH_OUTPUT_PATH, "MC", year, f"MC.{sample}_{year}")
            for variation in sample_variations:
                print("  >> " + variation)
                fpath = base_path + '_' + variation + ".root"
                self.set_sample_var_fields(year, sample, fpath, variation)

    def split_samples_by_flavour(self, sample):
        """
        We only treat the DY+Jets process split into light and heavy
        flavour. This function returns a list of suffixes corresponding
        to the light and heavy flavour part of the process.
        No suffix means inclusive treatment.
        """
        if "DYJets" in sample:
            return ["_bjet", "_ljet"]
        if "WJets" in sample:
            return ["_bjet", "_ljet"]
        else:
            return [""]

    def load_flavour_selection_mask(self, f, jet_flav):
        hadflav = f["AnalysisTree/slimmedGenJets.m_hadronFlavour"].array(library="np")
        gpt = f["AnalysisTree/slimmedGenJets.m_pt"].array(library="np")
        sel_bjet = np.array([np.any((i[j > 20] == 4) | (i[j > 20] == 5)) for i, j in zip(hadflav, gpt)])

        if jet_flav == "_bjet":
            sel_gen = sel_bjet
        elif jet_flav == "_ljet":
            sel_gen = ~sel_bjet
        else:
            sel_gen = None

        return sel_gen

    def set_sample_var_fields(self, year, sample, fpath, variation):
        fpath = fpath.replace("JES_","JES")
        fpath = fpath.replace("JER_","JER")
        with uproot.open(fpath) as f:
            for jet_flav in self.split_samples_by_flavour(sample):
                sel_gen = self.load_flavour_selection_mask(f, jet_flav)

                for x in self.variables_to_load:
                    if "ak4chs_Total_nominal" in x[0]:
                        continue
                    key, branch = self._tree_to_np_array(f, x)
                    if sel_gen is not None:
                        branch = branch[sel_gen]
                    pa_table = pa.table({"foo": branch})
                    key = key.replace("/",".")
                    pq.write_table(pa_table, f"{CACHE}/mc_{year}_{sample}{jet_flav}_{variation}_{key}.parquet")

    @property
    def _varied_weights_to_load(self):
        return [x + xvar for x, xvars in self.config.variations.items() for xvar in xvars]

    def _get_variable_key(self, x):
        if isinstance(x, list):
            return f"{x[0]}_{x[1]}"
        else:
            return x

    def load_trees(self, sample):

        for year in self.ul_years:
            print(f"Loading {year} {sample}")

            mc_path = os.path.join(self.UHH_OUTPUT_PATH, "MC", year, f"MC.{sample}_{year}.root")

            with uproot.open(mc_path) as f:
                for jet_flav in self.split_samples_by_flavour(sample):
                    sel_gen = self.load_flavour_selection_mask(f, jet_flav)

                    for x in self.variables_to_load:
                        key, branch = self._tree_to_np_array(f, x)
                        if sel_gen is not None:
                            branch = branch[sel_gen]
                        pa_table = pa.table({"foo": branch})
                        key = key.replace("/",".")
                        pq.write_table(pa_table, f"{CACHE}/mc_{year}_{sample}{jet_flav}_nominal_{key}.parquet")

                    # Variations
                    for branch in itertools.chain(self._varied_weights_to_load):
                        if(("pdf" in branch) and (("up" in branch) or ("down" in branch))):
                            continue
                        branch_name = re.sub(r"_mu[rf]{1}_", "_murmuf_", branch)
                        branch_name = re.sub(r"_[ab]{1}$", "", branch_name)
                        x = f["AnalysisTree/" + branch_name].array(library="np")
                        if sel_gen is not None:
                            x = x[sel_gen]
                        if "trigger" in branch:  # TODO: What's going on here?? What is this doing?
                            x = np.where(x != 0, x, 1)
                        pa_table = pa.table({"foo": x})
                        #branch = branch.replace("/",".")
                        pq.write_table(pa_table, f"{CACHE}/mc_{year}_{sample}{jet_flav}_variation_{branch}.parquet")

    def merge_16_pre_post(self):

        def concat_variation(sample, variation, x):
            print(f"{CACHE}/mc_UL16pre/postVFP_{sample}_{variation}_{x}.parquet")
            pq_table_pre = pq.read_table(f"{CACHE}/mc_UL16preVFP_{sample}_{variation}_{x}.parquet")
            pq_table_post = pq.read_table(f"{CACHE}/mc_UL16postVFP_{sample}_{variation}_{x}.parquet")
            return np.concatenate(
                (pq_table_pre["foo"].to_numpy(),
                 pq_table_post["foo"].to_numpy())
            )

        def concat_nominal(sample, variable):
            if sample == "data":
                mc_or_data = "data"
                postfix = variable
            else:
                mc_or_data = "mc"
                postfix = f"{sample}_nominal_{variable}"
            print(f"{CACHE}/{mc_or_data}_UL16preVFP_{postfix}.parquet")
            pq_table_pre = pq.read_table(f"{CACHE}/{mc_or_data}_UL16preVFP_{postfix}.parquet")
            pq_table_post = pq.read_table(f"{CACHE}/{mc_or_data}_UL16postVFP_{postfix}.parquet")
            return np.concatenate(
                (pq_table_pre["foo"].to_numpy(),
                 pq_table_post["foo"].to_numpy())
            )

        _samples = self.samples
        if _samples == ["DYJets"]:
            _samples = ["DYJets_bjet", "DYJets_ljet"]
        if _samples == ["WJets"]:
            _samples = ["WJets_bjet", "WJets_ljet"]
        print(_samples)
        for sample in _samples:
            # Concatenate Nominal Trees
            for variable in self.variables_to_load:
                variable = self._get_variable_key(variable)
                pa_table = pa.table({"foo": concat_nominal(sample, variable)})
                if sample == "data":
                    pq.write_table(pa_table, f"{CACHE}/data_UL16_{variable}.parquet")
                else:
                    pq.write_table(pa_table, f"{CACHE}/mc_UL16_{sample}_nominal_{variable}.parquet")

            if sample == "data":
                continue

            # Concatenate "nominal sample" Variations
            for variation in itertools.chain(self._varied_weights_to_load):
                if(("pdf" in variation) and (("up" in variation) or ("down" in variation))):
                    continue
                pa_table = pa.table({"foo": concat_variation(sample, "variation", variation)})
                pq.write_table(pa_table, f"{CACHE}/mc_UL16_{sample}_variation_{variation}.parquet")

            # Concatenate sample Variations
            for sample_variation in self._get_relevant_sample_variations(sample):
                for x in self.variables_to_load:
                    key = self._get_variable_key(x)
                    if "ak4chs_" in key:
                        continue
                    pa_table = pa.table({"foo": concat_variation(sample, sample_variation, key)})
                    key = key.replace("/",".")
                    pq.write_table(pa_table, f"{CACHE}/mc_UL16_{sample}_{sample_variation}_{key}.parquet")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", default="", help="Optional name of a"
                        "single signal (corresponding to one mA/mH mass point) for which"
                        "to run combine factory.")
    args = parser.parse_args()
    Parquetifier(args.sample)