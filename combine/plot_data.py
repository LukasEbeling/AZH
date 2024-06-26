#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import os
import sys
import argparse
import uproot
import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
import mplhep as hep
from math import sqrt

from combine import Combine
from plot_utils import PlotMeta, CMSSW_BASE, REGIONS, BACKGROUNDS

sys.path.append("../factory")
from utils import TEMPLATES

plt.style.use(hep.style.CMS)

BKG_LABEL_MAP = {
        "DYJets_bjet": r"Z/$\gamma$ + HF",
        "DYJets_ljet": r"Z/$\gamma$ + LF",
        "WJets_bjet": r"W + HF",
        "WJets_ljet": r"W + LF",
        "SingleTop": r"single $t$",
        "QCD": r"QCD",
        "VV": r"diboson",
        "TTW": r"$t \bar t W$",
        "TTZ": r"$t \bar t Z$",
        "TT": r"$t \bar t$",
}

class PlotMetaDataMC(PlotMeta):

    OBS_XLABEL_MAP = {
        "MET": r"$p_T^{miss}$ [GeV]",
        "MTA": r"$m_{T,A}$ [GeV]",
        "MH": r"$m_{H}$ [GeV]",
        "Jet1Phi": r"leading jet $\phi$",
        "Jet1Eta": r"leading jet $\eta$",
        "Jet1Pt": r"leading jet $p_T$ [GeV]",
        "Jet2Phi": r"subleading jet $\phi$",
        "Jet2Eta": r"subleading jet $\eta$",
        "Jet2Pt": r"subleading jet $p_T$ [GeV]",
        "HT": r"$H_T$ [GeV]",
        "2DEllipses": r"ellipsis",
        "num_b": r"N_b",
        "num_j": r"N_{jets}",
        "num_l": r"N_{lep}",
        "b_pt": r"highest b-score $p_T$",
        "b_eta": r"highest b-score $\eta$",
        "b_phi": r"highest b-score $\phi$",
    }



@dataclass
class Fitter():
 
    signal: str
    obs: str
    channel: str
    region: str
    fit: str
    unblind: bool
    hists: dict

    def run_fits(self):
        basename = f"AZH_{self.signal}_{self.obs}_{self.channel}_{self.region}"
        card = f"{TEMPLATES}/UL17/{basename}.dat"
        #workspace = f"tmp/{basename}.root"
        workspace = f"{basename}.root"
        #shapes = f"tmp/shapes.{basename}.root"
        shapes = f"shapes.{basename}.root"
        
        combine = Combine(initilize=True)
        combine.create_workspace(card,workspace)
        combine.fit_blind(card,workspace)
        combine.save_shapes(card,workspace,shapes)


    def load(self):
        basename = f"AZH_{self.signal}_{self.obs}_{self.channel}_{self.region}"
        shapes = f"tmp/shapes.{basename}.root"

        self.hists = {"data": None, "bkgs": {}}
        with uproot.open(shapes) as f:
            self.hists["data"] = f[f"{self.region}_{self.fit}/data_obs"]
            self.hists["total_bkg"] = f[f"{self.region}_{self.fit}/TotalBkg"]
            self.hists["AtoZH"] = f[f"{self.region}_{self.fit}/AtoZH"].to_numpy()
            for bkg in BACKGROUNDS:
                try:
                    self.hists["bkgs"][bkg] = f[
                        f"{self.region}_{self.fit}/{bkg}"
                    ].to_numpy()
                except KeyError:
                    thist = self.hists["total_bkg"].to_numpy()
                    self.hists["bkgs"][bkg] = np.histogram(-1 * thist[0] - 1, thist[1])


    def plot(self):
        plot_meta = PlotMetaDataMC()

        fig, ax = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(10, 10),
            gridspec_kw={"height_ratios": (4, 1)},
            sharex=True,
        )

        # Styling with mplhep
        hep.cms.label(
            ax=ax[0],
            llabel=plot_meta.LLABEL,
            data=True,
            lumi=plot_meta.lumi("UL17"),
            year=plot_meta.year("UL17"),
        )

        # Histogram
        bins = self.hists["data"].to_numpy()[1]
        bin_widths = np.diff(bins)

        error_band_args = {
            "edges": bins,
            "facecolor": "none",
            "linewidth": 0,
            "alpha": 0.9,
            "color": "black",
            "hatch": "///",
        }

        mc = self.hists["total_bkg"].to_numpy()[0]
        mc_error = self.hists["total_bkg"].errors() 
        # Set data to total bkg in SR if run in blinded mode
        if ("SR" not in self.region) or self.unblind:
            data = self.hists["data"].to_numpy()[0]
            data_error = self.hists["data"].errors()
        else:
            data = mc
            data_error = [sqrt(bin_entry) for bin_entry in data] #correct?
            data = data * (-100)

        hep.histplot(
            [x[0]/bin_widths for x in self.hists["bkgs"].values()],
            bins=bins,
            histtype="fill",
            stack=True,
            #label=list(self.hists["bkgs"].keys()),
            label = [BKG_LABEL_MAP[bkg] for bkg in self.hists["bkgs"].keys()],
            color=plot_meta.colors(self.hists["bkgs"].keys()),
            ax=ax[0],
        )
        ax[0].stairs(
            (mc_error + mc)/bin_widths,
            baseline=(mc - mc_error)/bin_widths,
            label=f"{self.fit} Uncertainty",
            **error_band_args,
        )
        if "SR" not in self.region:
            hep.histplot(
                data/bin_widths,
                bins=bins,
                yerr=data_error/bin_widths,
                histtype="errorbar",
                markersize=13,
                color="k",
                label="data",
                ax=ax[0],
            )
        mA = self.signal.split("_")[0]
        mH = self.signal.split("_")[1]
        if "SR" in self.region:
            hep.histplot(
                self.hists["AtoZH"][0]*5/bin_widths, #rescale to AZH->ttvv at 1pb
                bins=bins,
                histtype="step",
                yerr=False,
                label=rf"$m_{{A(H)}} = {mA} ({mH})$ GeV",
                linewidth=2.5,
                alpha=0.85,
                color="tab:red",
                ax=ax[0],
            )

        # Ratio
        pulls = False
        if pulls:
            sigma_tot = np.sqrt(
                mc_error**2 + data_error**2
            )  # TODO: Figure out how MC stats enteres
            plls = (data - mc) / sigma_tot
            ax[1].plot(
                bins[:-1] + (bins[1] - bins[0]) / 2,
                plls,
                color="k",
                marker="o",
                linestyle="None",
            )
            ax[1].axhline(0, color="grey", linestyle="solid")
            ax[1].axhline(1, color="grey", alpha=0.7, linestyle="dashed")
            ax[1].axhline(-1, color="grey", alpha=0.7, linestyle="dashed")
            ax[1].axhline(2, color="grey", alpha=0.4, linestyle="dotted")
            ax[1].axhline(-2, color="grey", alpha=0.4, linestyle="dotted")
            ax[1].set_ylabel(r"$\frac{Data-MC}{\sigma_{tot}}$")
            ax[1].set_ylim([-2.8, 2.8])
            ax[1].yaxis.set_major_locator(MaxNLocator(integer=True))
        else:
            hep.histplot(
                data / mc,
                yerr=data_error / mc,
                histtype="errorbar",
                bins=bins,
                stack=False,
                color="k",
                ax=ax[1],
            )
            ax[1].stairs(
                1 + mc_error / mc, baseline=1 - mc_error / mc, **error_band_args
            )
            ax[1].axhline(1, color="grey", linestyle="--")
            ax[1].set_ylabel("data/pred",fontsize=25)
            ax[1].set_ylim([0.5, 1.5])

        # Plot Settings
        #title = rf"$\it{{{self.region}}}$  $\it{{{self.channel}}}$ $\it{{Channel}}$"
        title = REGIONS[self.region]
        ax[0].set_yscale("log")
        handles, labels = ax[0].get_legend_handles_labels()
        legend = ax[0].legend(
            reversed(handles),
            reversed(labels),
            loc="upper right",
            ncol=3,
            title=title,
            fontsize=18,
            title_fontsize=18,
            frameon=True,
        )
        #legend.get_frame().set_facecolor('white')
        if self.obs == "2DEllipses": ax[0].set_ylabel("events/bin",fontsize=25)
        else: ax[0].set_ylabel("events/GeV",fontsize=25)
        ax[0].set_ylim(ymin=1e-2, ymax=np.max(mc/bin_widths) * 900)
        if self.obs in ["beta","bphi","Jet1Eta","Jet2Eta","Jet1Phi","Jet2Phi","2DEllipses"]: 
            ax[0].set_ylim(ymin=1e-2, ymax=np.max(mc/bin_widths) * 90000)
        #ax[0].set_ylim(ymin=0, ymax=np.max(mc) + 100)
        max_bin = bins[-1] if bins[-1] < 10000 else bins[-2] #crop overflow bin
        ax[0].set_xlim(bins[0],max_bin) 
        ax[1].set_xlim(bins[0],max_bin)
        ax[1].set_xlabel(plot_meta.xlabel(self.obs))

        ax[0].grid()
        ax[1].grid()
        ax[0].tick_params(axis='y', labelsize=20)
        ax[1].tick_params(axis='y', labelsize=20)
        ax[1].tick_params(axis='x', labelsize=20, pad=10)

        # Save Plot
        fpath_out = os.path.join(
            CMSSW_BASE,
            f"src/UHH2/AZH/combine/plots/",
        )
        os.makedirs(fpath_out, exist_ok=True)
        fname = (
            self.channel
            + "_"
            + self.signal
            + "_"
            + self.obs
            + "_"
            + self.region
            + "_"
            + self.fit
        )
        if pulls:
            fname += "_pulls"
        plt.savefig(fpath_out + f"{fname}.png",bbox_inches="tight")
        #plt.savefig(fpath_out + f"{fname}.pdf",bbox_inches="tight")
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='eg ./plot_data 1000_400 MET --prefit')
    parser.add_argument('signal', type=str)
    parser.add_argument('obs', type=str)
    parser.add_argument('--prefit', action='store_true')

    args = parser.parse_args()
    signal = args.signal
    observable = args.obs
    fitting = "prefit" if args.prefit else "postfit"

    for region in REGIONS.keys():
        fitter = Fitter(
            signal=signal,
            #obs="MET",
            obs=observable,
            channel="inv",
            region=region,
            unblind=False,
            fit=fitting,
            hists={},
        )

        #fitter.run_fits()
        fitter.load()
        fitter.plot()
