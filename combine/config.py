#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import os
import yaml

from utils import CMSSW_BASE


class Borg():

    _shared_state = {}

    def __init__(self):
        self.__dict__ = self._shared_state


class Configurator(Borg):

    def __init__(self, signal: str = ""):
        # Make Configuration Singleton
        Borg.__init__(self)
        if self._shared_state:
            return

        self.signal = signal
        with open(os.path.join(CMSSW_BASE, "src/UHH2/AZH/combine/config.yaml"), 'r') as f:
            config_dict = yaml.safe_load(f)
        self.config_dict = config_dict

    @property
    def years(self):
        return self.config_dict["years"]

    @property
    def signals(self):
        if self.signal:
            return ["AZH_" + self.signal]

        with open(os.path.join(CMSSW_BASE, "src/UHH2/AZH/config/signals.txt"), 'r') as f:
            return [f"AZH_{line}".removesuffix("\n").replace("MA-","").replace("MH-","") for line in f]

    @property
    def backgrounds(self):
        return self.config_dict["backgrounds"]

    @property
    def samples(self):
        return self.signals + self.backgrounds

    @property
    def svars(self):
        return self.config_dict["svars"]

    @property
    def branches(self):
        return self.config_dict["branches"]

    @property
    def region_branch(self):
        return self.config_dict["region_branch"]

    @property
    def variations(self):
        return self.config_dict["variations"]

    @property
    def sample_variations(self):
        return self.config_dict["sample_variations"]

    @property
    def sample_var_whitelist(self):
        sample_var_whitelist = self.config_dict["sample_var_whitelist"]
        for sample_variation in sample_var_whitelist:
            if sample_var_whitelist[sample_variation] == "ALL PROCESSES":
                sample_var_whitelist[sample_variation] = self.samples + ["DYJets", "AtoZH"]
        return sample_var_whitelist

    @property
    def angle_cut_on(self):
        return self.config_dict["angle_cut"]["active"]

    @property
    def angle_cut_signal_efficiency(self):
        return self.config_dict["angle_cut"]["signal_efficiency"]

    @property
    def dnn_cut_on(self):
        return self.config_dict["dnn_cut"]["active"]

    @property
    def dnn_cut_signal_efficiency(self):
        return self.config_dict["dnn_cut"]["signal_efficiency"]


if __name__ == "__main__":
    c = Configurator()
    print(c.variations)
    print("Signals: ", c.signals)
    print("Backgrounds: ", c.backgrounds)
    print("Years: ", c.years)
    print(c.svars)
    print(c.angle_cut_on)
    print(c.angle_cut_signal_efficiency)
    print(c.dnn_cut_on)
    print(c.dnn_cut_signal_efficiency)
    print(c.sample_var_whitelist)
    print(' ')
    print(c.sample_variations.items())