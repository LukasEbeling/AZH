#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

from os.path import isfile, join, exists
from os import environ, mkdir
from datetime import datetime
from itertools import product
from progress.bar import IncrementalBar
import uproot
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import sys

#plt.style.use(hep.style.CMS)

CMSSW_BASE = environ.get('CMSSW_BASE')
FILES_PATH = join(CMSSW_BASE, "src/UHH2/AZH/data/output_02_reconstruction/MC/UL17/")
OUTPUT_PATH = join(CMSSW_BASE,"src/UHH2/AZH/plots")
BACKGROUNDS = ["QCD","SingleTop", "TT","TTZ", "WJets", "TTW", "VV","ZJets"]
SIGNALS = ["1000_400","600_400","700_450","750_400","750_650","800_400","1000_400","1000_850"]
CUTS = ["base", "met", "phi", "btag", "veto", "weight"]
BRANCH_MAP = {
    "base": "CutFlow_Baseline",
    "met": "CutFlow_MET>170",
    "phi": "CutFlow_deltaphi",
    "btag": "CutFlow_TwoB",
    "veto": "CutFlow_LeptonVeto",
    "weight": "CutFlow_Weights"
}
LABEL_MAP = {
    "base": "base",
    "met": "met$\geq$170",
    "phi": "$\Delta \phi\geq$0.5",
    "btag": "#b$\geq$2",
    "veto": "veto",
    "weight": "weight$\leq$10"
}

class DataLoader():
    def __init__(self):
        self.signal = {}
        self.background = {}
        self.load_bkgs()
        self.load_sign()

    def load_bkgs(self):
        for bkg in BACKGROUNDS: 
            path = join(CMSSW_BASE,FILES_PATH,f"MC.{bkg}_UL17.root")
            self.background[bkg] = self.load(path)
    
    def load_sign(self):
        for sig in SIGNALS: 
            path = join(CMSSW_BASE,FILES_PATH,f"MC.INV_{sig}_UL17.root")
            self.signal[sig] = self.load(path)

    def load(self, path):
        counts = []
        with uproot.open(path) as file:
            for cut in CUTS:
                branch = BRANCH_MAP[cut]
                events = file[branch]["N_Events"].to_numpy()
                counts.append(float(events[0][0]))
        return np.array(counts)

def plot(axis,ratio,name):
    edges = np.linspace(1,len(ratio)+1,len(ratio)+1)
    
    axis.hist(edges[:-1], weights=ratio, bins=edges, histtype='step',
        label=name)

    x_coordinates = np.arange(0, len(ratio)) + 1.5
    x_labels = [LABEL_MAP[cut] for cut in CUTS]
    axis.xaxis.set_major_locator(plt.FixedLocator(x_coordinates))
    axis.xaxis.set_major_formatter(plt.FixedFormatter(x_labels))
    axis.set_yscale('log') 
    axis.legend(loc='lower left')
    plt.title("Cutflow")


def make_figure(data,legend,title):
    fig, ax = plt.subplots(figsize=(7,3))
    doc = open(OUTPUT_PATH+f"/yields_reco_{title}.txt", "w")

    for name in legend:
        counts = data[name]
        plot(ax,counts/counts[0],name)
        doc.write(name+" "+str(counts[-1])+"\n")
        bar.next()

    fig.savefig(join(OUTPUT_PATH,f"cutflow_reco_{title}"))
    plt.close(fig)
    doc.close()


if __name__ == "__main__":
    data_loader = DataLoader()
    bar = IncrementalBar(f"Making Plots", max=len(BACKGROUNDS)+len(SIGNALS))
    data = data_loader.background
    make_figure(data,BACKGROUNDS,"bkg")
    data = data_loader.signal
    make_figure(data,SIGNALS,"sig")
    bar.finish()
