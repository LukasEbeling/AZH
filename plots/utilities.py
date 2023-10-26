#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import json
import os
import re
import numpy as np

ANALYSIS = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH"


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
        return get_inclusive(mA,mH) * br_Z_vv


class Limit:

    def __init__(self, channel="inv", var="MET", region="all", year="UL17"):
        self.channel = channel
        self.var = var
        self.region = region
        self.year = year

    def load(self, mA, mH, pct = "50.0%"):
        fname = f"AZH_{mA}_{mH}_{self.var}_{self.channel}_{self.region}.log"
        with open(os.path.join(ANALYSIS, "combine", self.year, fname), "r") as f:
            limit_line = [x for x in f if pct in x][0]

        limit_value = float(re.findall(r"\d+.\d+", limit_line.split("<")[1])[0])
        return limit_value


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


