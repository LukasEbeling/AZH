#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import os

import matplotlib.pyplot as plt


def CMSSW_BASE():
    return "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28"


class PlotMeta:

    YEAR_LUMI_MAP = {
        "UL18": 59.8,
        "UL17": 41.5,
        "UL16": 36.3,
        "ULCombined": 138
    }

    UL_YEAR_MAP = {
        "UL18": 2018,
        "UL17": 2017,
        "UL16": 2016,
    }

    CMAP = plt.cm.get_cmap("Set3")

    LLABEL = "Private work"

    def year(self, ul_year: str):
        return self.UL_YEAR_MAP[ul_year]

    def lumi(self, ul_year: str):
        return self.YEAR_LUMI_MAP[ul_year]

    def xlabel(self, obs: str):
        return self.OBS_XLABEL_MAP[obs]

    def colors(self, processes: list):
        id_map = {
            "AtoZH": 8,
            #"DYJets": 0,
            "ZJets": 0,
            "TT": 1,
            "TTZ": 2,
            "QCD": 3,
            "VV": 4,
            "SingleTop": 5,
            "TTW": 6,
            "WJets": 7,
        }

        return [self.CMAP.colors[id_map[p]] for p in processes]


class map_region_int_to_str:
    single_region_map = {
        "all": "All Control Regions",
        0: "SignalRegion",
        5: "SRFiveJets",
        51: "CRFiveJets1BTagged",
        50: "CRFiveJets0BTagged",
        11: "CRZMassSidebands",
        110: "CRZMassSidebands0BTagged",
        111: "CRZMassSidebands1BTagged",
        512: "CRFiveJetsZMassSidebands",
        510: "CRFiveJetsZMassSidebands0BTagged",
        511: "CRFiveJetsZMassSidebands1BTagged",
        12: "CRSameSignLeptons",
        121: "CRSameSignLeptons1BTagged",
        120: "CRSameSignLeptons0BTagged",
        13: "CR0BTaggedJets",
        14: "CR1BTaggedJets",
        15: "CRDiffLeptonFlavours"
    }

    def __getitem__(self, i):
        if isinstance(i, str) or isinstance(i, int):
            return self.single_region_map[i]
        if isinstance(i, list):
            return '\n+'.join([self.single_region_map[x] for x in i])


map_year_lumi = {
    'UL18': 59.83,
    'UL17': 41.48,
    'UL16preVFP': 19.5,
    'UL16postVFP': 16.8,
}


map_ul_year = {
    'UL18': 2018,
    'UL17': 2017,
    'UL16preVFP': 2016,
    'UL16postVFP': 2016,
}

map_channel_int_to_str = {
    11 * 11: "Electron",
    13 * 13: "Muon",
    11 * 13: "EMu"
}


def is_valid_combination(channel, region, obs):
    if channel == 11 * 11 and "muon" in obs:
        return False
    if channel == 13 * 13 and "electron" in obs:
        return False
    if channel == 11 * 13 and region != 15:
        return False
    if channel != 11 * 13 and region == 15:
        return False
    return True


def errorband(
    H_band,  # value array if error band is symmetric / tuple of upper and lower value arrays
    bins,  # corresponding to bin edges of the hist values supplied in H_band
    label=None,
    edges=True,
    ax=None,
    **kwargs,
):
    # ax check
    if ax is None:
        ax = plt.gca()
    else:
        if not isinstance(ax, plt.Axes):
            raise ValueError("ax must be a matplotlib Axes object")

    # Construction of upper and lower boundries
    y_upper = H_band[0]
    y_lower = H_band[1]

    error_band_args = {
        "edges": bins, "facecolor": "none", "linewidth": 0,
        "alpha": .9, "color": "black", "hatch": "///"
    }
    ax.stairs(y_upper, baseline=y_lower, label=label, **{**error_band_args, **kwargs})


if __name__ == "__main__":
    pass