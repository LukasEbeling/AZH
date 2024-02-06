#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import os
import json
import re
import numpy as np

import matplotlib.pyplot as plt


CMSSW_BASE = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28"
ANALYSIS = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH"

BACKGROUNDS = [
    "WJets_bjet", 
    "WJets_ljet", 
    "DYJets_bjet", 
    "DYJets_ljet", 
    "VV", 
    "QCD", 
    "SingleTop",
    "TTW", 
    "TTZ", 
    "TT",
]

REGIONS = {
    "SR_2B_6J" : r"SR 0l-2b-6j",
    "SR_2B_5J" : r"SR 0l-2b-5j", 
    "SR_1B_6J" : r"SR 0l-1b-6j",
    "SR_1B_5J" : r"SR 0l-1b-5j",
    "IR_0B_5J" : r"CR 0l-0b-5j",
    "IR_0B_6J" : r"CR 0l-0b-6j",
    "LR_2B_5J" : r"CR 1l-2b-5j",  
    "LR_2B_6J" : r"CR 1l-2b-6j",
    "LR_1B_5J" : r"CR 1l-1b-5j",
    "LR_1B_6J" : r"CR 1l-1b-6j",
    "LR_0B_5J" : r"CR 1l-0b-5j",
    "LR_0B_6J" : r"CR 1l-0b-6j",
    #"SR_DNN": r"SR node",
    #"TT_DNN": r"TT node",
    #"QCD_DNN": r"QCD node",
    #"DY_DNN": r"DYJets node",
    #"WJ_DNN": r"WJets node",
}

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

    #CMAP = plt.cm.get_cmap("Set1")
    CMAP = plt.cm.get_cmap("Set3")

    LLABEL = "Private work"

    def year(self, ul_year: str):
        return self.UL_YEAR_MAP[ul_year]

    def lumi(self, ul_year: str):
        return self.YEAR_LUMI_MAP[ul_year]

    def xlabel(self, obs: str):
        return self.OBS_XLABEL_MAP[obs]

    def colors(self, processes: list):
        id_map_alt = {
            "AtoZH": 5,
            "DYJets": 0,
            "DYJets_ljet": 0,
            "DYJets_bjet": 0,
            "TT": 8,
            "TTZ": 2,
            "QCD": 6,
            "VV": 1,
            "SingleTop": 7,
            "TTW": 4,
            "WJets": 3,
            "WJets_ljet": 3,
            "WJets_bjet": 3,
        }

        id_map = {
            "AtoZH": 0,
            "TT": 1,
            "TTZ": 13,
            "TTW": 9,
            "SingleTop": 19,
            "QCD": 6,
            "VV": 4,
            "DYJets": 0,
            "DYJets_ljet": 7,
            "DYJets_bjet": 18,
            "WJets": 0,
            "WJets_ljet": 15,
            "WJets_bjet": 2,
        }
        CMAP1 = plt.cm.get_cmap("Set1")
        CMAP2 = plt.cm.get_cmap("Set3")
        color = CMAP1.colors + CMAP2.colors
        return [color[id_map[p]] for p in processes]


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

class Theory: 

    def __init__(self,tanb = 1):
        tanb_str = str(float(tanb))

        with open(f"theory_predictions/xsec_{tanb_str}.json", "r") as f:
            self.xsecs = json.load(f)

        with open(f"theory_predictions/brs_{tanb_str}.json", "r") as f:
            self.brs = json.load(f)


    def get_inclusive(self,mA, mH): 
        br_A_ZH = self.brs[f"{float(mA)},{float(mH)}"]["a_zh"]
        br_H_tt = self.brs[f"{float(mA)},{float(mH)}"]["h_tt"]
        xsec = self.xsecs[f"{float(mA)}"]
        return xsec * br_A_ZH * br_H_tt

    def get_invisible(self,mA,mH):
        br_Z_vv = 0.2
        return self.get_inclusive(mA,mH) * br_Z_vv

class Limit:

    def __init__(self, channel="inv", var="MET", region="all", year="UL17"):
        self.channel = channel
        self.var = var
        self.region = region
        self.year = year

    def load(self, mA, mH, pct = "50.0%"):
        file = f"tmp/AZH_{mA}_{mH}_{self.var}_{self.channel}_{self.region}.log"
        with open(file, "r") as f:
            limit_line = [x for x in f if pct in x][0]

        limit_value = float(re.findall(r"\d+.\d+", limit_line.split("<")[1])[0])
        return limit_value

class BestLimit:

    def __init__(self, channel="inv", region="all", year="UL17"):
        self.channel = channel
        self.region = region
        self.year = year

    def load(self, mA, mH, pct = "50.0%"):
        channel = self.channel
        region = self.region
        year = self.year

        median_met = Limit(channel,"MET",region,year).load(mA,mH,"50.0%")
        median_mh = Limit(channel,"MH",region,year).load(mA,mH,"50.0%")
        median_mt = Limit(channel,"MTA",region,year).load(mA,mH,"50.0%")
        median_2d = Limit(channel,"2DEllipses",region,year).load(mA,mH,"50.0%")
        index = np.argmin((median_met,median_mh,median_mt,median_2d),axis=0)
        limit_met = Limit(channel,"MET",region,year).load(mA,mH,pct)
        limit_mh = Limit(channel,"MH",region,year).load(mA,mH,pct)
        limit_mt = Limit(channel,"MTA",region,year).load(mA,mH,pct)
        limit_2d = Limit(channel,"2DEllipses",region,year).load(mA,mH,pct)
        return (limit_met,limit_mh,limit_mt,limit_2d)[index]

def load_masses():
    with open(
        os.path.join(ANALYSIS, "config/signals.txt"), "r"
    ) as f:
        mass_points = np.array(
            [
                (int(re.findall(r"\d+", x)[0]), int(re.findall(r"\d+", x)[1]))
                for x in f
            ]
        )
    #mass_points = np.array([(MA,MH) for MA,MH in mass_points if MA<1001 and MH<1001])
    return mass_points