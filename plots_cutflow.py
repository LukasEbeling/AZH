#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

from os.path import isfile, join, exists
from os import environ, mkdir
from datetime import datetime
from itertools import product
import uproot
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import sys

#plt.style.use(hep.style.CMS)

CMSSW_BASE = environ.get('CMSSW_BASE')
OUTPUT_PATH = join(CMSSW_BASE, "src/UHH2/AZH/data/output_01_preselection/MC/UL17/")
BACKGROUNDS = ["QCD"]
SIGNALS = ["1000_400","600_400","700_450","750_400","750_650","800_400","1000_400","1000_850"]
CUTS = ["base","veto","50met","100met","150met","6j","phi","2b"]
BRANCH_MAP = {
    "base": "CutFlow_Baseline",
    "veto": "CutFlow_LeptonVeto",
    "50met": "CutFlow_MET>50",
    "100met": "CutFlow_MET>100",
    "150met": "CutFlow_MET>150",
    "6j": "CutFlow_SixJets",
    "phi": "CutFlow_deltaphi",
    "2b": "CutFlow_TwoB",
}

class DataLoader():
    def __init__(self):
        self.signal = {}
        self.background = {}
        self.load_bkgs()
        self.load_sign()

    def load_bkgs(self):
        for bkg in BACKGROUNDS: 
            path = join(CMSSW_BASE,OUTPUT_PATH,f"MC.{bkg}_UL17.root")
            self.background[bkg] = self.load(path)
    
    def load_sign(self):
        for sig in SIGNALS: 
            path = join(CMSSW_BASE,OUTPUT_PATH,f"MC.INV_{sig}_UL17.root")
            self.signal[sig] = self.load(path)

    def load(self, path):
        cutflow = []
        with uproot.open(path) as file:
            for cut in CUTS:
                branch = BRANCH_MAP[cut]
                events = file[branch]["N_Events"].to_numpy()
                cutflow.append(events[0][0])
        return cutflow/cutflow[0]

def plot(axis,counts,name):
    edges = np.linspace(1,len(counts)+1,len(counts)+1)
    
    axis.hist(edges[:-1], weights=counts, bins=edges, histtype='step',
        label=name)

    x_coordinates = np.arange(0, len(counts)) + 1.5
    axis.xaxis.set_major_locator(plt.FixedLocator(x_coordinates))
    axis.xaxis.set_major_formatter(plt.FixedFormatter(CUTS))
    axis.set_yscale('log') 
    axis.legend()
    plt.title("Cutflow")


if __name__ == "__main__":
    data_loader = DataLoader()
    folder = join(CMSSW_BASE,"src/UHH2/AZH/plots")
    f = open(folder+"/yields.txt", "w")

    fig, ax = plt.subplots(figsize=(10,5))
    for bkg in BACKGROUNDS:
        counts = data_loader.background[bkg]
        plot(ax,counts,bkg)
        f.write(bkg+" "+str(counts[-1])+"\n")

    fig.savefig(join(folder,"cutflow_bkg"))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10,5))
    for sig in SIGNALS:
        counts = data_loader.signal[sig]
        plot(ax,counts,sig)
        f.write(sig+" "+str(counts[-1])+"\n")

    fig.savefig(join(folder,"cutflow_sig"))
    plt.close(fig)
    f.close()




