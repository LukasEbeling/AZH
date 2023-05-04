#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse
import itertools
import os
import re

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import uproot

from config import Configurator
import utils


CMSSW_BASE = os.environ.get('CMSSW_BASE')


class Parquetifier():

    UHH_OUTPUT_PATH = os.path.join(CMSSW_BASE, "src/UHH2/AZH/data/output_02_reconstruction/")

    def __init__(self, sample=None):
        self.config = Configurator()
        self.ul_years = utils.split_ul16(self.config.years)
        self.samples = [re.sub(r"^MA-", "INV_", sample)]
        if not sample:
            self.samples = self.config.samples
        self.nominal_trees = {year: {} for year in self.ul_years}
        self.sample_vars = {year: {} for year in self.ul_years}

        self.vars_to_load = self.config.svars + self.config.branches

        os.makedirs("cache", exist_ok=True)

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
                for x in self.vars_to_load:
                    key, branch = self._tree_to_np_array(f, x)
                    pa_table = pa.table({"foo": branch})
                    pq.write_table(pa_table, f"cache/data_{year}_{key}.parquet")

    def load_sample_vars(self, sample):
        if self.config.sample_variations is None: return #remove line later
        for year in self.ul_years:
            print(f"Loading Jet Variations {year} {sample}")
            relevant_vars = {
                x: y for x, y in
                self.config.sample_variations.items() if
                sample in self.config.sample_var_whitelist[x]
            }
            sample_variations = [x + xvar for x, xvars in relevant_vars.items() for xvar in xvars]

            base_path = os.path.join(self.UHH_OUTPUT_PATH, "MC", year, f"MC.{sample}_{year}")
            for variation in sample_variations:
                print("  >> " + variation)
                fpath = base_path + '_' + variation + ".root"
                self.set_sample_var_fields(year, sample, fpath, variation)

    def set_sample_var_fields(self, year, sample, fpath, variation):
        with uproot.open(fpath) as f:
            for x in self.vars_to_load:
                key, branch = self._tree_to_np_array(f, x)
                pa_table = pa.table({"foo": branch})
                pq.write_table(pa_table, f"cache/mc_{year}_{sample}_{variation}_{key}.parquet")

    def load_trees(self, sample):

        for year in self.ul_years:
            print(f"Loading {year} {sample}")

            mc_path = os.path.join(self.UHH_OUTPUT_PATH, "MC", year, f"MC.{sample}_{year}.root")
            with uproot.open(mc_path) as f:
                for x in self.vars_to_load:
                    key, branch = self._tree_to_np_array(f, x)
                    pa_table = pa.table({"foo": branch})
                    pq.write_table(pa_table, f"cache/mc_{year}_{sample}_nominal_{key}.parquet")

                    
    def load_trees_var(self, sample):

        for year in self.ul_years:
            print(f"Loading {year} {sample}")

            mc_path = os.path.join(self.UHH_OUTPUT_PATH, "MC", year, f"MC.{sample}_{year}.root")
            with uproot.open(mc_path) as f:
                for x in self.vars_to_load:
                    key, branch = self._tree_to_np_array(f, x)
                    pa_table = pa.table({"foo": branch})
                    pq.write_table(pa_table, f"cache/mc_{year}_{sample}_nominal_{key}.parquet")

                # Variations
                weights_to_load = [x + xvar for x, xvars in
                                   self.config.variations.items() for xvar in xvars]

                for branch in itertools.chain(weights_to_load):
                    if(("pdf" in branch) and (("up" in branch) or ("down" in branch))):
                        continue
                    branch_name = re.sub(r"_mu[rf]{1}_", "_murmuf_", branch)
                    branch_name = re.sub(r"_[ab]{1}$", "", branch_name)
                    x = f["AnalysisTree/" + branch_name].array(library="np")
                    if "trigger" in branch:  # TODO: What's going on here?? What is this doing?
                        x = np.where(x != 0, x, 1)
                    pa_table = pa.table({"foo": x})
                    pq.write_table(pa_table, f"cache/mc_{year}_{sample}_variation_{branch}.parquet")

    def merge_16_pre_post(self):  # TODO: Fix this function

        def concat_sample_var(sample, variation, x):
            return np.concatenate(
                (self.sample_vars["UL16preVFP"][sample][variation][x],
                 self.sample_vars["UL16postVFP"][sample][variation][x])
            )

        def concat_nominal(sample, x):  # TODO: load arrays from pq files here and concatenate those
            return np.concatenate(
                (self.nominal_trees["UL16preVFP"][sample][x],
                 self.nominal_trees["UL16postVFP"][sample][x])
            )

        # Concatenate Nominal Trees
        for sample in self.nominal_trees["UL16preVFP"].keys():
            self.nominal_trees["UL16"][sample] = {}
            for x in list(self.nominal_trees["UL16preVFP"][sample].keys()):  # TODO: Adapt this loos
                pa_table = pa.table({"foo": concat_nominal(sample, x)})
                pq.write_table(pa_table, f"cache/mc_UL16_{sample}_nominal_{x}.parquet")

        for sample in self.samples:
            # Concatenate Jet Variations
            for variation in self.sample_vars["UL16preVFP"][sample].keys():
                for x in self.sample_vars["UL16preVFP"][sample][variation].keys():
                    print(sample, variation, x)
                    pa_table = pa.table({"foo": concat_sample_var(sample, variation, x)})
                    pq.write_table(pa_table, f"cache/mc_UL16_{sample}_{variation}_{x}.parquet")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", default="", help="Optional name of a"
                        "single signal (corresponding to one mA/mH mass point) for which"
                        "to run combine factory.")
    args = parser.parse_args()
    Parquetifier(args.sample)
