#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import itertools
import os
import sys

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

sys.path.append("../factory/")
from generate_datacards import SHAPES_NP
from generate_datacards import RATE_NP
from config import Configurator
from utils import CMSSW_BASE, TEMPLATES

plt.style.use(hep.style.CMS)

OUTPUT_PATH = os.path.join(CMSSW_BASE, "src/UHH2/AZH/combine/nuisances")

config = Configurator()
SIGNALS = ["AZH_1000_400"]
OBSERVABLES = ["MET"]
EXCLUDE_BKGS = []
CHANNELS = ["inv"]
REGIONS = ["SR_2B_6J"]
YEARS = ['UL17']
YEAR_LUMI_MAP = {
    "UL18": 59.83,
    "UL17": 41.48,
    "UL16": 36.3,
}
UL_YEAR_MAP = {
    "UL18": 2018,
    "UL17": 2017,
    "UL16": 2016,
}
OBS_XLABEL_MAP = {
    "2DEllipses": r"$p_{T,mis} \times m_H$",
    "MET": r"$p_{T,mis}$",
    "DeltaM": r"$\Delta m = m_A - m_H$",
    "ZPT": r"$p_{T}(Z)$",
    "Elec1Pt": r"$p_{T}(e_1)$",
    "Muon1Pt": r"$p_{T}(\mu_1)$",
    "Jet1Pt": r"$p_{T}(j_1)$",
    "ZMass": r"$m_Z$",
}


def load(_year, _signal, _obs, _ch, _reg):
    """
    Loads needed TH1s from the combine input files.
    """
    fname = f"{_signal}_{_obs}_{_ch}_{_reg}.root"
    hists = {"nominals": {}}
    with uproot.open(os.path.join(TEMPLATES, _year, fname)) as f:
        # Load nominal
        available_processes = []

        # Load all processes that are available
        for process in config.samples + ["AtoZH"]:
            try:
                if process in EXCLUDE_BKGS:
                    raise uproot.exceptions.KeyInFileError(process)
                hists["nominals"][process] = f[f"{_reg}/{process}"]
                available_processes.append(process)
            except uproot.exceptions.KeyInFileError:
                pass

        #SHAPES_NP.update(RATE_NP)
        for shape_np, np_processes in SHAPES_NP.items():
            shape_np = shape_np.replace("YEAR",_year)
            hists[shape_np] = {}
            for process in np_processes:
                if not process in available_processes:
                    print(process)
                    continue
                hists[shape_np][process] = {}
                hists[shape_np][process]["up"] = f[f"{_reg}/{process}_{shape_np}Up"]
                hists[shape_np][process]["down"] = f[f"{_reg}/{process}_{shape_np}Down"]
    return hists

def save_plot(path, file):
    os.makedirs(path, exist_ok=True)
    plt.savefig(os.path.join(path,file+".png"))

def plotUncertaintyRatios(hists, _year, _signal, _obs, _ch, _reg, _process, _np):
    """
    Really ugly plotting function.
    """
    no_color = "k"
    up_color = "tab:red"
    dn_color = "tab:blue"
    var_linestyle = "dashed"
    lw = 2

    fig, ax = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(12, 10),
        gridspec_kw={"height_ratios": (4, 1)},
        sharex=True,
    )

    hep.cms.label(
        ax=ax[0],
        llabel="Work in progress",
        data=True,
        lumi=YEAR_LUMI_MAP[_year],
        year=UL_YEAR_MAP[_year],
    )

    # Histogram
    bins = hists[_np][_process]["up"].to_numpy()[1]
    error_band_args = {
        "edges": bins,
        "facecolor": "none",
        "linewidth": 0,
        "alpha": 0.74,
    }
    nominal = hists["nominals"][_process].to_numpy()[0]
    var_up = hists[_np][_process]["up"].to_numpy()[0]
    var_down = hists[_np][_process]["down"].to_numpy()[0]
    hep.histplot(
        nominal,
        yerr=0,
        bins=bins,
        linewidth=lw,
        color=no_color,
        histtype="step",
        label="nominal",
        ax=ax[0],
    )
    hep.histplot(
        var_up,
        bins=bins,
        linestyle=var_linestyle,
        linewidth=lw,
        color=up_color,
        histtype="step",
        label="up",
        ax=ax[0],
    )
    hep.histplot(
        var_down,
        bins=bins,
        linestyle=var_linestyle,
        linewidth=lw,
        color=dn_color,
        histtype="step",
        label="down",
        ax=ax[0],
    )

    # Ratio
    hep.histplot(
        var_up / nominal,
        histtype="step",
        bins=hists["nominals"][_process].to_numpy()[1],
        stack=False,
        linestyle=var_linestyle,
        color=up_color,
        ax=ax[1],
    )
    hep.histplot(
        var_down / nominal,
        histtype="step",
        bins=hists["nominals"][_process].to_numpy()[1],
        stack=False,
        linestyle=var_linestyle,
        color=dn_color,
        ax=ax[1],
    )
    ax[1].axhline(1, color=no_color)
    ax[1].set_ylabel(r"$\frac{Variation}{Nominal}$")
    ax[1].set_ylim([0.6, 1.4])

    # Plot Settings
    title = f"{_process} {_np}\n{_reg}, {_year}, \n{_ch} Channel {OBS_XLABEL_MAP[_obs]}"
    ax[0].set_yscale("log")
    ax[0].legend(ncol=2, title=title, fontsize=18, title_fontsize=18, frameon=False)
    ax[0].set_ylabel("Events")
    ax[0].set_ylim(ymin=1e-1)
    ax[1].set_xlabel(OBS_XLABEL_MAP[_obs])

    # Save Plot
    fname = (
        "UncertaintyRatio"
        + "_"
        + _np
        + "_"
        + _ch
        + "_"
        + _reg
        + "_"
        + _obs
        + "_"
        + _year
        + "_"
        + _signal
        + "_"
        + _process
    )

    fpath_np = os.path.join(OUTPUT_PATH, f"{_np}/")
    fpath_np = fpath_np.replace("_"+process,"")
    save_plot(fpath_np, fname)
    plt.close()


if __name__ == "__main__":

    for year, signal, obs, ch, reg in itertools.product(
        YEARS, SIGNALS, OBSERVABLES, CHANNELS, REGIONS
    ):
        hists = load(year, signal, obs, ch, reg)
        for nparam in list(set(hists.keys()) - set(["nominals"])):
            available_processes = list(
                set.intersection(
                    set(hists[nparam].keys()), set(hists["nominals"].keys())
                )
            )
            for process in available_processes:
                if process in EXCLUDE_BKGS:
                    continue
                #print(year, signal, obs, ch, reg, process, nparam)
                plotUncertaintyRatios(
                    hists, year, signal, obs, ch, reg, process, nparam
                )