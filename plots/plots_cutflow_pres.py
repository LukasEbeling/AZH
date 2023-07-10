#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

from os.path import isfile, join, exists
from os import environ, mkdir
from datetime import datetime
from itertools import product
from operator import add
from progress.bar import IncrementalBar
import uproot
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import sys
import glob

#plt.style.use(hep.style.CMS)

CMSSW_BASE = environ.get('CMSSW_BASE')
FILES_PATH = join(CMSSW_BASE, "src/UHH2/AZH/data/output_01_preselection/MC/UL17/")
OUTPUT_PATH = join(CMSSW_BASE,"src/UHH2/AZH/plots")
BACKGROUNDS = []
SIGNALS = ["1000_400","600_400","700_450","750_400","750_650","800_400","1000_400","1000_850"]
CUTS = ["base", "met", "jets"]
BRANCH_MAP = {
    "base": "CutFlow_Baseline",
    "met": "CutFlow_MET>50",
    "jets": "CutFlow_SixJets",
}
LABEL_MAP = {
    "base": "base",
    "met": "met$\geq$50",
    "jets": "#jets$\geq$6",
}

class DataLoader():
    def __init__(self):
        self.signal = {}
        self.background = {}
        self.load_bkgs()
        self.load_sign()

        self.background["QCD"] = self.load_all("QCD")
        self.background["ST"] = self.load_all("ST")
        self.background["TTZ"] = self.load_all("TTZ")
        self.background["TTW"] = self.load_all("TTW")
        self.background["TT"] = self.load_all("TTo")
        self.background["WJets"] = self.load_all("WJets")
        self.background["ZJets"] = [x+y for x,y in zip(self.load_all("ZJets"),self.load_all("DY"))]
        self.background["VV"] = [x+y+z for x,y,z in zip(self.load_all("WW"),self.load_all("ZZ"),self.load_all("WZ"))]

        #doc = open(OUTPUT_PATH+f"/yields_preselection.txt", "w")
        #for sample in self.signal:
        #    doc.write("\""+sample+"\": "+str(list(self.signal[sample]))+",\n")
        #for sample in self.background:
        #    doc.write("\""+sample+"\": "+str(list(self.background[sample]))+",\n")
        #doc.close()

    
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

    def load_all(self,name):
        paths = glob.glob(join(FILES_PATH, f"*{name}*"))
        #if not isfile(path): continue 
        counts = np.zeros(len(CUTS))
        for path in paths:
            counts = list(map(add,counts,self.load(path))) 
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


def make_figure(data,title):
    fig, ax = plt.subplots(figsize=(7,3))
    #doc = open(OUTPUT_PATH+f"/yields_{title}.txt", "w")

    names = data.keys()
    for name in names:
        counts = data[name]
        plot(ax,counts/counts[0],name)
        #doc.write(name+" "+str(counts[-1])+"\n")
        bar.next()

    fig.savefig(join(OUTPUT_PATH,f"cutflow_{title}"))
    plt.close(fig)
    #doc.close()


if __name__ == "__main__":
    data_loader = DataLoader()
    bkg_data = data_loader.background
    sig_data = data_loader.signal
    bar = IncrementalBar(f"Making Plots", max=len(bkg_data.keys())+len(sig_data.keys()))
    make_figure(bkg_data,"bkg")
    make_figure(sig_data,"sig")
    bar.finish()
