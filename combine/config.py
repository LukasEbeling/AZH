#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import os
import yaml


CMSSW_BASE = os.environ.get("CMSSW_BASE")


class Borg():

    _shared_state = {}

    def __init__(self):
        self.__dict__ = self._shared_state


class Configurator(Borg):

    def __init__(self):
        Borg.__init__(self)

        if self._shared_state:
            return

        with open(os.path.join(CMSSW_BASE, "src/UHH2/AZH/combine/config.yaml"), 'r') as f:
            config_dict = yaml.safe_load(f)
        with open(os.path.join(CMSSW_BASE, "src/UHH2/AZH/config/signals.txt"), 'r') as f:
            self.signals = [f"INV_{line}".removesuffix("\n").replace("MA-","").replace("MH-","") for line in f]
        print(self.signals)

        self.years = config_dict["years"]
        self.backgrounds = config_dict["backgrounds"]
        self.samples = self.signals + self.backgrounds
        self.svars = config_dict["svars"]
        self.branches = config_dict["branches"]
        #self.region_branch = config_dict["region_branch"]
        self.variations = config_dict["variations"]
        #self.sample_variations = config_dict["sample_variations"]
        #self.sample_var_whitelist = config_dict["sample_var_whitelist"]
        #for jes_variation in ["Total", "Absolute", "Absolute_YEAR", "BBEC1", "BBEC1_YEAR", "EC2", "EC2_YEAR", "FlavorQCD", "HF", "HF_YEAR", "RelativeBal", "RelativeSample_YEAR"]:
        #    self.sample_var_whitelist[jes_variation] = config_dict["sample_var_whitelist"][jes_variation] + self.signals
        #self.sample_var_whitelist["jer"] = config_dict["sample_var_whitelist"]["jer"] + self.signals
        #self.angle_cut_on = config_dict["angle_cut"]["active"]
        #self.angle_cut_signal_efficiency = config_dict["angle_cut"]["signal_efficiency"]
        #self.dnn_cut_on = config_dict["dnn_cut"]["active"]
        #self.dnn_cut_signal_efficiency = config_dict["dnn_cut"]["signal_efficiency"]


if __name__ == "__main__":
    c = Configurator()
    print("Signals: ", c.signals)
    print("Backgrounds: ", c.backgrounds)
    print("Years: ", c.years)
    print(c.svars)

