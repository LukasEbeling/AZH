#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse

#from elliptical_hist_plotter import EllipticalHistPlotter
from ntuple_loader import NTupleLoader
import utils
from utils import REGION_ID_MAP
from sample_handler import SampleSet
from config import Configurator


class UHH2ToCombineFactory():
    """
    Creates root files in the structure required for
    combine fitting procedures with histograms binned
    accoring to the specifications of our analysis.
    This class is the high level factory that triggers
    the creation of all desired sample sets.
    """

    def __init__(self, signal: str = ""):
        config = Configurator()
        if signal:
            signal = "AZH_" + signal
            config.signals = [signal]
            config.samples = [signal] + config.backgrounds
            self.loader = NTupleLoader(signal)

        svars = list(map(utils.collection_key, config.svars))
        # create a dictionary containing a keys and SampleSets
        channel = "inv"
        region = "SignalRegion"

        self.sample_sets = {
            f"{year}_{channel}_{region}_{sgnl}_{svar}": SampleSet(
                {
                    "signal": sgnl,
                    "channel": channel,
                    "region": region,
                    "year": year,
                    "svar": svar
                }
            )
            for sgnl in config.signals
            for year in config.years
            for svar in svars
            if utils.is_valid_set(channel, region, svar)
        }

    def run_factory(self):
        # First iteration without 2D CRs
        for set_name, sample_set in self.sample_sets.items():
            in_signal_region = sample_set.set_params["region"] == "SignalRegion"
            if (in_signal_region):
                sample_set.set_sample_bins()
                sample_set.create_hists_and_save_file()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--signal", default="", help="Optional name of a"
                        "single signal (corresponding to one mA/mH mass point) for which"
                        "to run combine factory.")
    args = parser.parse_args()
    file_factory = UHH2ToCombineFactory(args.signal)
    file_factory.run_factory()

