import itertools
import os
import re

import pyarrow.parquet as pq

from config import Configurator
from utils import CMSSW_BASE,CACHE


class BorgNTupleLoader():

    _shared_state = {}

    def __init__(self):
        self.__dict__ = self._shared_state


class NTupleLoader(BorgNTupleLoader):

    def __init__(self, signal: str = ""):
        # Makse NTupleLoader Singleton
        BorgNTupleLoader.__init__(self)
        if self._shared_state:
            return

        self.config = Configurator(signal)
        self.nominal_trees = {year: {} for year in self.config.years}
        self.sample_vars = {year: {} for year in self.config.years}

        self.vars_to_load = self.config.svars + self.config.branches

        for s in self.config.samples:
            self.load_trees(s)
            self.load_sample_vars(s)
        self.load_data()

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
                pq_table = pq.read_table(f"{CACHE}/data_{year}_{key}.parquet")
                self.nominal_trees[year]["data"][key] = pq_table["foo"].to_numpy()

    def load_sample_vars(self, sample):

        for year in self.sample_vars:
            # Uncomment the following lines to produce jet related cards
            # if not self.sample_vars[year]:
            #     continue
            print(f"Loading Jet Variations {year} {sample}")
            relevant_vars = {
                x: y for x, y in
                self.config.sample_variations.items() if
                sample in self.config.sample_var_whitelist[x]
            }
            sample_variations = [x + xvar for x, xvars in relevant_vars.items() for xvar in xvars]
            self.sample_vars[year][sample] = {vari: {} for vari in sample_variations}

            for variation in sample_variations:
                print("  >> " + variation)
                self.set_sample_var_fields(year, sample, variation)

    def set_sample_var_fields(self, year, sample, variation):
        for x in self.vars_to_load:
            key = self._tree_to_np_array(x)
            pq_table = pq.read_table(f"{CACHE}/mc_{year}_{sample}_{variation}_{key}.parquet")
            self.sample_vars[year][sample][variation][key] = pq_table["foo"].to_numpy()

    def load_trees(self, sample):

        for year in self.nominal_trees:
            print(f"Loading {year} {sample}")

            self.nominal_trees[year][sample] = {}
            for x in self.vars_to_load:
                key = self._tree_to_np_array(x)
                # self.nominal_trees[year][sample][key] = branch
                pq_table = pq.read_table(f"{CACHE}/mc_{year}_{sample}_nominal_{key}.parquet")
                self.nominal_trees[year][sample][key] = pq_table["foo"].to_numpy()

            # Variations
            weights_to_load = [x + xvar for x, xvars in
                               self.config.variations.items() for xvar in xvars]

            for branch in itertools.chain(weights_to_load):
                if(("pdf" in branch) and (("up" in branch) or ("down" in branch))):
                    continue
                branch_name = re.sub(r"_mu[rf]{1}_", "_murmuf_", branch)
                branch_name = re.sub(r"_[ab]{1}$", "", branch_name)
                pq_table = pq.read_table(f"{CACHE}/mc_{year}_{sample}_variation_{branch}.parquet")
                self.nominal_trees[year][sample][branch] = pq_table["foo"].to_numpy()