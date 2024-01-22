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


    def __init__(self, signal: str = ""):
        config = Configurator(signal)
        if signal:
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
            for channel in ["inv"]
            for region in REGION_ID_MAP.keys()
            for year in config.years
            for svar in svars
            for svar in svars + ["ellipses"]
            #if utils.is_valid_set(channel, region, svar)
        }

    def _run_sr_ellipses_in_crs(self):
        # First iteration without 2D CRs
        for set_name, sample_set in self.sample_sets.items():
            in_signal_region = ("SR" in sample_set.set_params["region"])
            no_2d_binning = ("ellipse" not in sample_set.set_params["svar"])
            if (no_2d_binning or in_signal_region):
                sample_set.set_sample_bins()
                sample_set.create_hists_and_save_file()
        # Second iteration only 2D CRs fetching the ellipses from
        # SR sample sets that are already finished
        for set_name, sample_set in self.sample_sets.items():
            in_control_region = ("SR" not in sample_set.set_params["region"])
            elliptical_binning = ("ellipse" in sample_set.set_params["svar"])
            if (elliptical_binning and in_control_region):                
                control_region = sample_set.set_params["region"]
                signal_region = "SR_2B_5J" if "5J" in control_region else "SR_2B_6J"
                #signal_region = control_region.replace("LR","SR").replace("IR","SR")                
                sr_set_name = set_name.replace(control_region, signal_region) 
                ellipses = self.sample_sets[sr_set_name].set_binning
                sample_set.set_sample_bins(ellipses)
                sample_set.create_hists_and_save_file()

    def _run_simple(self):
        for set_name, sample_set in self.sample_sets.items():
            sample_set.set_sample_bins()
            sample_set.create_hists_and_save_file()
    
    def run_factory(self):
        self._run_sr_ellipses_in_crs()
        #or self._run_simpe()
 
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
    signal = args.signal
    signal = signal.replace("MA-","")
    signal = signal.replace("MH-","")

    file_factory = UHH2ToCombineFactory(signal)
    file_factory.run_factory()
    file_factory.plot_elliptical_binnings()