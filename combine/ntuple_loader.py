import itertools
import os
import re

import pyarrow.parquet as pq
import numpy as np

from config import Configurator
import utils


CMSSW_BASE = os.environ.get('CMSSW_BASE')


class BorgNTupleLoader():

    _shared_state = {}

    def __init__(self):
        self.__dict__ = self._shared_state


class NTupleLoader(BorgNTupleLoader):

    UHH_OUTPUT_PATH = os.path.join(CMSSW_BASE, "src/UHH2/AZH/data/output_02_reconstruction/")

    def __init__(self, signal: str = ""):
        BorgNTupleLoader.__init__(self)
        if self._shared_state:
            return

        self.config = Configurator()
        ul_years = utils.split_ul16(self.config.years)
        if signal:
            self.samples = [signal] + self.config.backgrounds
        else:
            self.samples = self.config.samples
        self.nominal_trees = {year: {} for year in ul_years}
        self.sample_vars = {year: {} for year in ul_years}

        self.vars_to_load = self.config.svars + self.config.branches

        for s in self.samples:
            self.load_trees(s)
            self.load_sample_vars(s)
        self.load_data()

        if "UL16" in self.config.years:
            self.merge_16_pre_post()

    def _tree_to_np_array(self, x):
        # This is the case in which the observable
        # is part of a collection, e.g. jet or lepton pt
        if isinstance(x, list):
            key = f"{x[0]}_{x[1]}"
        else:
            key = x
        return key

    def load_data(self):

        for year in self.nominal_trees:
            print(f"Loading {year} DATA")
            self.nominal_trees[year]["data"] = {}
            for x in self.vars_to_load:
                key = self._tree_to_np_array(x)
                pq_table = pq.read_table(f"cache/data_{year}_{key}.parquet")
                self.nominal_trees[year]["data"][key] = pq_table["foo"].to_numpy()

    def load_sample_vars(self, sample):
        if self.config.sample_variations is None: return #remove line later

        for year in self.sample_vars:
            print(f"Loading Jet Variations {year} {sample}")
            relevant_vars = {
                x: y for x, y in
                self.config.sample_variations.items() if
                sample in self.config.sample_var_whitelist[x]
            }
            sample_variations = [x + xvar for x, xvars in relevant_vars.items() for xvar in xvars]
            self.sample_vars[year][sample] = {vari: {} for vari in sample_variations}

            base_path = os.path.join(self.UHH_OUTPUT_PATH, "MC", year, f"MC.{sample}_{year}")
            for variation in sample_variations:
                print("  >> " + variation)
                fpath = base_path + '_' + variation + ".root"
                self.set_sample_var_fields(year, sample, fpath, variation)

    def set_sample_var_fields(self, year, sample, fpath, variation):
        for x in self.vars_to_load:
            key = self._tree_to_np_array(x)
            pq_table = pq.read_table(f"cache/mc_{year}_{sample}_{variation}_{key}.parquet")
            self.sample_vars[year][sample][variation][key] = pq_table["foo"].to_numpy()

    def load_trees(self, sample):

        for year in self.nominal_trees:
            print(f"Loading {year} {sample}")

            self.nominal_trees[year][sample] = {}
            for x in self.vars_to_load:
                key = self._tree_to_np_array(x)
                ## self.nominal_trees[year][sample][key] = branch
                pq_table = pq.read_table(f"cache/mc_{year}_{sample}_nominal_{key}.parquet")
                self.nominal_trees[year][sample][key] = pq_table["foo"].to_numpy()

            # Variations
            weights_to_load = [x + xvar for x, xvars in
                               self.config.variations.items() for xvar in xvars]

            for branch in itertools.chain(weights_to_load):
                if(("pdf" in branch) and (("up" in branch) or ("down" in branch))):
                    continue
                branch_name = re.sub(r"_mu[rf]{1}_", "_murmuf_", branch)
                branch_name = re.sub(r"_[ab]{1}$", "", branch_name)
                pq_table = pq.read_table(f"cache/mc_{year}_{sample}_variation_{branch}.parquet")
                self.nominal_trees[year][sample][branch] = pq_table["foo"].to_numpy()

    def merge_16_pre_post(self):

        def concat_sample_var(sample, variation, x):
            return np.concatenate(
                (self.sample_vars["UL16preVFP"][sample][variation][x],
                 self.sample_vars["UL16postVFP"][sample][variation][x])
            )

        def concat_nominal(sample, x):
            return np.concatenate(
                (self.nominal_trees["UL16preVFP"][sample][x],
                 self.nominal_trees["UL16postVFP"][sample][x])
            )

        self.nominal_trees["UL16"] = {}
        self.sample_vars["UL16"] = {}

        # Concatenate Nominal Trees
        for sample in self.nominal_trees["UL16preVFP"].keys():
            self.nominal_trees["UL16"][sample] = {}
            for x in list(self.nominal_trees["UL16preVFP"][sample].keys()):
                pq_table = pq.read_table(f"cache/mc_UL16_{sample}_nominal_{x}.parquet")
                self.nominal_trees["UL16"][sample][x] = pq_table["foo"].to_numpy()

        for sample in self.samples:
            # Concatenate Jet Variations
            self.sample_vars["UL16"][sample] = {}

            for variation in self.sample_vars["UL16preVFP"][sample].keys():
                self.sample_vars["UL16"][sample][variation] = {}
                for x in self.sample_vars["UL16preVFP"][sample][variation].keys():
                    print(sample, variation, x)
                    pq_table = pq.read_table(f"cache/mc_UL16_{sample}_{variation}_{x}.parquet")
                    self.sample_vars["UL16"][sample][variation][x] = pq_table["foo"].to_numpy()
