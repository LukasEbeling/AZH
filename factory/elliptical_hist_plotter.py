import os
import sys

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

sys.path.append("../combine/")
import plot_utils  # noqa


CMSSW_BASE = os.environ.get("CMSSW_BASE")
OUTPUT_PATH = os.path.join(CMSSW_BASE, "src/UHH2/AZH/factory/plots/")
os.makedirs(OUTPUT_PATH, exist_ok=True)


class EllipticalHistPlotter():

    ell_colors = plt.cm.tab20b((4 / 3 * np.arange(20 * 3 / 4)).astype(int))
    sgnl_color = "tab:red"

    def __init__(self, sample_sets):
        self.sample_sets = sample_sets
        #for s in sample_sets:
        #    for ellipses in s.set_binning:
        #        print(ellipses.n_std)

    def plot(self):
        for sample_set in self.sample_sets:
            self._plot_plane_with_bins(sample_set)
            #self._plot_twod_unrolled_bins(sample_set)

    def _plot_skeleton(self, ul_year):
        fig, ax = plt.subplots(dpi=300)
        lumi = plot_utils.PlotMeta.YEAR_LUMI_MAP[ul_year]
        year = plot_utils.PlotMeta.UL_YEAR_MAP[ul_year]
        hep.cms.label("Private work", lumi=lumi, year=year, data=True, loc=0)
        plt.style.use(hep.style.CMS)
        return fig, ax

    def _plot_plane_with_bins(self, sample_set):
        sgnl_raw = sample_set.samples['signal']._get_raw_twod_inputs()
        bkg_met = []
        bkg_dm = []
        for sname, s_bkg in sample_set.samples.items():
            if sname in ["signal", "data"]:
                continue
            met, dm = s_bkg._get_raw_twod_inputs()
            bkg_met.append(met)
            bkg_dm.append(dm)

        fig = plt.figure(figsize=(11, 11), dpi=300)

        # 1D Plots
        left, width = 0.15, 0.6
        bottom, height = 0.15, 0.6
        spacing = 0.00

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom + height + spacing, width, 0.15]
        rect_histy = [left + width + spacing, bottom, 0.15, height]

        ax = fig.add_axes(rect_scatter)
        ax_histx = fig.add_axes(rect_histx, sharex=ax, frameon=True)
        ax_histy = fig.add_axes(rect_histy, sharey=ax, frameon=True)
        tick_args = {
            "reset": True,
            "labelbottom": False,
            "labelleft": False
        }
        ax_histx.tick_params(**tick_args)
        ax_histy.tick_params(**tick_args)
        ax.tick_params(direction='in')

        ax.hexbin(sgnl_raw[0], sgnl_raw[1], gridsize=35, cmap='Blues',
                  vmin=.1)
        ax.set_facecolor(plt.get_cmap("Blues")(0))
        for ell, color in zip(sample_set.set_binning, self.ell_colors):
            ax.add_patch(ell.get_plt_patch(color))

        # 2D Plots
        binwidth = 25
        xymax = np.max(np.abs(sgnl_raw[0]))
        lim = (int(xymax / binwidth) + 1) * binwidth

        bins = np.arange(-lim, lim + binwidth, binwidth)
        ax_histx.hist(sgnl_raw[0], bins=bins)
        ax_histy.hist(sgnl_raw[1], bins=2 * bins, orientation='horizontal')

        ax.set_xlim(0)
        ax.set_ylim(50)
        ax.set_xlabel(r"$p_T^{miss}$ [GeV]", size=25)
        ax.set_ylabel(r"$m_H$ [GeV]", size=25)
        ax.legend(
            [ell.get_plt_patch(c) for ell, c in zip(sample_set.set_binning, self.ell_colors)],
            [rf"{ell.n_std}$\sigma$" for ell in sample_set.set_binning],
            fontsize=20
        )
        ax.tick_params(axis='both', which='major', labelsize=18)
        ul_year = sample_set.set_params["year"]
        fname_out = os.path.join(
            OUTPUT_PATH,
            "ellipses_"
            + sample_set.samples["signal"].sample
            + "_"
            + sample_set.set_params["region"]
            + f"_{ul_year}"
        )
        #plt.savefig(fname_out + ".png")
        plt.savefig(fname_out + ".pdf")
        plt.close()

    def _plot_twod_unrolled_bins(self, sample_set):
        sgnl, sgnl_raw, bins = sample_set.samples['signal'].svar_hist
        bkgs = []
        bkg_labels = []
        for sname, bkg_sample in sample_set.samples.items():
            if sname in ["signal", "data"]:
                continue
            bkgs.append(bkg_sample.svar_hist[0])
            bkg_labels.append(bkg_sample.sample)

        plot_meta = plot_utils.PlotMeta()

        fig, ax = self._plot_skeleton(ul_year=sample_set.set_params["year"])
        ax2 = ax.twinx()
        ax2.hist(bins[:-1], bins=bins, weights=sgnl,
                 label=sample_set.samples["signal"].sample, histtype='step', linewidth=3,
                 color=self.sgnl_color, linestyle=(0, (5, 5)))
        ax2.set_ylabel('# Signal Events', color=self.sgnl_color)
        ax2.tick_params(axis='y', labelcolor=self.sgnl_color)
        ax.set_ylabel("# Background Events", size=30)
        ax.set_xlabel("Standard Deviations", size=30)
        #bins = [i for i in range(len(bins))] if "overflow" in bins else bins
        ax.hist([[x - .01 for x in bins[1:]] for _ in bkgs], bins=bins, weights=bkgs,
                label=bkg_labels, stacked=True,
                color=plot_meta.colors(bkg_labels))
        ax.legend(loc="lower left")
        fig.tight_layout()
        ul_year = sample_set.set_params["year"]
        fname_out = os.path.join(
            OUTPUT_PATH, "unrolled_bins_"
            + sample_set.samples["signal"].sample
            + "_"
            + sample_set.set_params["region"]
            + f"_{ul_year}"
        )
        plt.savefig(fname_out + ".png")
        plt.savefig(fname_out + ".pdf")
        plt.close()