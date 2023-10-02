#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
from collections import OrderedDict as od
from itertools import product
from math import pi
import os 

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pandas as pd
from progress.bar import IncrementalBar
import uproot


plt.style.use(hep.style.CMS)

#Config
ANALYSIS = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/"
RECO_PATH = ANALYSIS + "data/output_02_reconstruction/"
CACHE = ANALYSIS + "combine/cache/"
BACKGROUNDS = ["VV", "TTW", "TTZ", "DYJets_ljet", "DYJets_bjet", "WJets_ljet", "WJets_bjet", "QCD", "SingleTop","TT"]
MASSES = ["600_400","1000_400","2000_400"]
SIGNALS = [f"AZH_{mass}" for mass in MASSES]
CMAP = plt.cm.get_cmap("Set1")

REGIONS = {
    0: "SR_6J",
    1: "SR_5J",
    11: "IR_1B_5J",
    12: "IR_1B_6J",
    13: "IR_0B_5J",
    14: "IR_0B_6J",
    15: "LR_2B_5J",
    16: "LR_2B_6J",
    17: "LR_1B_5J",
    18: "LR_1B_6J",
    19: "LR_0B_5J",
    20: "LR_0B_6J",
}

OBSERVABLES = {
    "MET": np.linspace(170,800,10),
    #"delta_m": np.linspace(200,1500,13),
    #"delta_mt": np.linspace(200,1500,13),
    "m_H": np.linspace(200,1500,13),
    #"mt_A": np.linspace(500,1800,13),
    #"mt_H": np.linspace(200,1500,13),
    #"jetsAk4CHS/jetsAk4CHS.m_pt": np.linspace(80, 800, 21),
    #"jetsAk4CHS/jetsAk4CHS.m_phi": np.linspace(-pi, pi, 21),
    #"jetsAk4CHS/jetsAk4CHS.m_eta": np.linspace(-3, 3, 21),
    #"event_weight": np.linspace(0,100,21), 
}

COLOR = {
    "1000_400": "black",
    "DYJets": CMAP.colors[0],
    "DYJets_ljet": CMAP.colors[0],
    "DYJets_bjet": CMAP.colors[0],
    "TT": CMAP.colors[8],
    "TTZ": CMAP.colors[2],
    "QCD": CMAP.colors[6],
    "VV": CMAP.colors[1],
    "SingleTop": CMAP.colors[7],
    "TTW": CMAP.colors[4],
    "WJets": CMAP.colors[3],
    "WJets_ljet": CMAP.colors[3],
    "WJets_bjet": CMAP.colors[3],
}




class dataLoader():

    data: dir = {}

    def __init__(self):
        #self.load_bkgs()
        #self.load_sigs()
        self.load_samples(BACKGROUNDS)
        self.load_samples(SIGNALS)

    def load(self, path):
        o = {}
        with uproot.open(path) as f:
            for obs in OBSERVABLES:
                o[obs] = f[f'AnalysisTree/{obs}'].array(library="np")
            w = f[f"AnalysisTree/event_weight"].array(library="np")
            r = f[f"AnalysisTree/region"].array(library="np")
        return o, w, r

    def load_new(self,sample):
        o = {}
        for obs in OBSERVABLES: o[obs] = self.load_parquet(sample,obs)
        w = self.load_parquet(sample, "event_weight")
        r = self.load_parquet(sample, "region")
        return o,w,r

    def load_parquet(self, sample, observable):
        file_path = CACHE + f"mc_UL17_{sample}_nominal_{observable}.parquet"
        df = pd.read_parquet(file_path)
        return np.array(df["foo"])
    
    def load_samples(self,samples):
        bar = IncrementalBar("Loading", max=len(samples))
        for sample in samples:
            print("\n"+sample)
            #path = RECO_PATH + "MC/UL17/" + f"MC.{sample}_UL17.root"
            #o, w, r = self.load(path)
            o,w,r = self.load_new(sample)
            self.data[sample] = {
                "observables": o,
                "weights": w,
                "regions": r,
            }
            bar.next()
        bar.finish()

def filter_events(val,wgh,reg,region,observable):
    v = val[reg == region]
    w = wgh[reg == region]

    if not len(v):
        assert not len(w)
        return v, w

    if "jets" in observable: v = np.array([jet[0] for jet in v])

    for i in range(len(w)): 
        if w[i]>10 and region==0: w[i]=0

    assert (v.shape == w.shape)
    return v, w

def overflow(val,binning):
    return np.clip(val,binning[0],binning[-1])

def add_bin(binning):
    return np.append(binning,binning[-1]+(binning[-1]-binning[-2]))


def get_bkg(region,binning,observable):  #as list 
    binned = []
    errors = []
 
    for bkg in BACKGROUNDS:
        val = data_loader.data[bkg]["observables"][observable]
        wgh = data_loader.data[bkg]["weights"]
        reg = data_loader.data[bkg]["regions"]

        val,wgh = filter_events(val,wgh,reg,region,observable)

        #met = data_loader.data[bkg]["observables"]["MET"]
        #met = met[reg==region]
        #for i in range(len(wgh)):
        #    if met[i]<50: wgh[i]=0

        val = overflow(val,binning)

        num,bins = np.histogram(val, weights=wgh, bins=binning)
        err,bins = np.histogram(val, weights=wgh**2, bins=binning)

        binned.append(num)
        errors.append(np.sqrt(err))

    return binned, errors

def get_sig(region,signal,binning,observable): 

    val = data_loader.data[signal]["observables"][observable]
    wgh = data_loader.data[signal]["weights"]
    reg = data_loader.data[signal]["regions"]

    val,wgh = filter_events(val,wgh,reg,region,observable)  
    val = overflow(val,binning)      

    num,bins = np.histogram(val, weights=wgh, bins=binning)
    err,bins = np.histogram(val, weights=wgh**2, bins=binning)

    err = np.sqrt(err)
    return num, err

def get_xlabel(obs):
    x = {
        "MET": (r"missing $p_T$ [GeV]","MET"),
        "m_H": (r"m of H [GeV]","MH"),
        "mt_A": (r"$m_T$ of A [GeV]","MTA"),
        "mt_H": (r"$m_T$ of H [GeV]","MTH"),
        "delta_m": (r"$m_T(A)-m(H)$","DM"),
        "delta_mt": (r"$m_T(A)-m_T(H)$","DMT"),
        "jetsAk4CHS/jetsAk4CHS.m_pt": (r"$p_{T}(j_{1})$ [GeV]", "Jet1PT"),
        "jetsAk4CHS/jetsAk4CHS.m_phi": (r"$\phi(j_{1})$ [rad]", "Jet1Phi"),
        "jetsAk4CHS/jetsAk4CHS.m_eta": (r"$\eta(j_{1})$", "Jet1Eta"),
        "event_weight": (r"Event Weight","Weight"),
        "W_m": (r"best $W_m$ [GeV]","Wmass"),
        "HT": (r"$H_T$ [GeV]","HT"),
        "delta_phi": (r"min $\Delta\phi$ [rad]","DeltaPhi"),
        "num_leptons": ("number of leptons", "LEP"),
        }
    return x[obs.replace("Puppi", "CHS")]

def plot_samples(region, signal, observable, binning=np.linspace(50, 750, 70)):
    binning = add_bin(binning)
    binned_bkg, errors = get_bkg(region,binning,observable)
    binned_sig, error = get_sig(region,signal,binning,observable)

    fig, axes = plt.subplots(
        nrows=2, 
        ncols=1,
        figsize=(12,10),
        gridspec_kw={"height_ratios": (4, 1)},
        sharex=True
    )

    yields_bkg = [round(sum(bkg),1) for bkg in binned_bkg]
    yield_sig = round(sum(binned_sig),1)

    label_bkg = [f"{b} ({y})" for b,y in zip(BACKGROUNDS,yields_bkg)]
    label_sig = f"{signal} ({yield_sig})" 

    hep.cms.label(ax=axes[0],llabel='Work in progress',data=True, lumi=41.48, year=2017)

    colors = [COLOR[process] for process in BACKGROUNDS]

    hep.histplot(binned_bkg, histtype='fill', bins=binning, stack=True, label=label_bkg, ax=axes[0],color=colors)
    hep.histplot(binned_sig, yerr=error, histtype='step', bins=binning, color='k', label=label_sig, ax=axes[0])


    #axes[0].set_yscale('log')
    axes[0].legend(ncol=2, title=REGIONS[region], fontsize=18,title_fontsize=18, frameon=True)
    axes[0].set_ylabel("Events")

    axes[1].set_xlabel(get_xlabel(observable)[0])

    # Save Plot
    #folder = os.path.join(ANALYSIS,"plots",signal)
    folder = os.path.join(ANALYSIS,"plots",signal)
    if not os.path.exists(folder): os.makedirs(folder)

    image = REGIONS[region].replace(' ', '') + '_' + get_xlabel(observable)[1] + "_UL17"

    plt.savefig(os.path.join(folder,image))
    plt.close()

if __name__ == "__main__":
    data_loader = dataLoader()
    data = data_loader.data

    bar = IncrementalBar("Plotting", max=len(OBSERVABLES.keys())*len(REGIONS.keys())*(len(SIGNALS)))

    for obs, binning in OBSERVABLES.items():
        for region, signal in product(REGIONS.keys(), SIGNALS):
            if region == 5 and obs == "MET": binning = np.linspace(50,150,10)
            plot_samples(region, signal, obs, binning)
            bar.next()
    bar.finish()
