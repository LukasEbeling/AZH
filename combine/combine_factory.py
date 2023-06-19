#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse

from elliptical_hist_plotter import EllipticalHistPlotter
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

    CHANNELS_REGIONS = {
        "inv": ["SignalRegion","CR_1L","CR_0B","CR_lowmet"],
    }

    def __init__(self, signal: str = ""):
        config = Configurator()
        if signal:
            signal = "AZH_" + signal
            config.signals = [signal]
            config.samples = [signal] + config.backgrounds
            self.loader = NTupleLoader(signal)

        svars = list(map(utils.collection_key, config.svars))
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
            for channel in self.CHANNELS_REGIONS
            for region in self.CHANNELS_REGIONS[channel]
            for year in config.years
            for svar in svars
            #for svar in svars + ["ellipses"]
            #if utils.is_valid_set(channel, region, svar)
        }

    def run_factory(self):
        # First iteration without 2D CRs
        for set_name, sample_set in self.sample_sets.items():
            sample_set.set_sample_bins()
            sample_set.create_hists_and_save_file()
        #If binning exists already:
        #ellipses = self.sample_sets[set_name].set_binning
        #sample_set.set_sample_bins(ellipses)

    def plot_elliptical_binnings(self):
        samples = [x for x in self.sample_sets.values() if x.set_params["svar"] == "ellipses"]
        ell_hist_plotter = EllipticalHistPlotter(samples)
        ell_hist_plotter.plot()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--signal", default="", help="Optional name of a"
                        "single signal (corresponding to one mA/mH mass point) for which"
                        "to run combine factory.")
    args = parser.parse_args()
    file_factory = UHH2ToCombineFactory(args.signal)
    file_factory.run_factory()
    #file_factory.plot_elliptical_binnings()