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

ANALYSIS = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/"
RECO_PATH = os.path.join(ANALYSIS, "data/output_02_reconstruction/")
BACKGROUNDS = ["QCD","SingleTop", "TT","TTZ", "WJets", "TTW", "VV"]
SIGNALS = ["1000_400"] #["600_400","700_450","750_400","750_650","800_400","1000_400","1000_850"]
REGIONS = [0] #0 corresponds to signalregion
YEARS = ['UL17']
YEAR_LUMI_MAP = {
    'UL18': 59.83,
    'UL17': 41.48,
    'UL16preVFP': 19.5,
    'UL16postVFP': 16.8,
    }
UL_YEAR_MAP = {
    'UL18': 2018,
    'UL17': 2017,
    'UL16preVFP': 2016,
    'UL16postVFP': 2016,
    }


class DataLoader():
    
    def __init__(self, observables):
        self.observables = observables
        self.sig = {}
        self.bkg = {}

        self.bar = IncrementalBar("Loading Samples", max=(len(SIGNALS)+len(BACKGROUNDS))*len(YEARS)*len(REGIONS))
        self.load_bkgs()
        self.load_sign()
        self.bar.finish()


    def load(self, _path):
        o = {}
        with uproot.open(_path) as f:
            for _obs in self.observables:
                o[_obs] = f[f'AnalysisTree/{_obs}'].array(library="np")
            w = f[f"AnalysisTree/event_weight"].array(library="np")
            r = f[f"AnalysisTree/region"].array(library="np")
        self.bar.next()
        return o, w, r


    def load_bkgs(self):
        for y in YEARS:
            self.bkg[y] = {}
            for bkg in BACKGROUNDS:
                _path = os.path.join(RECO_PATH, "MC", y, f"MC.{bkg}_{y}.root")
                o, w, r = self.load(_path)
                self.bkg[y][bkg] = {
                    "observables": o,
                    "weights": w,
                    "regions": r,
                    }
    

    def load_sign(self):
        for y in YEARS:
            self.sig[y] = {}
            for sig in SIGNALS:
                _path = os.path.join(RECO_PATH, "MC", y, f"MC.INV_{sig}_{y}.root")
                o, w, r = self.load(_path)
                self.sig[y][sig] = {
                    "observables": o,
                    "weights": w,
                    "regions": r,
                    }


def filter_events(observable, weights, region, _region, _obs):
    sel = region==_region
    o = observable[sel]
    w = weights[sel]

    if not len(o):
        assert not len(w)
        return o, w

    if "muon" in _obs or "electron" in _obs:
        # Only consider leading lepton information
        o = np.array([j[0] for j in o])
    if "jets" in _obs:
        # Only consider leading jet information
        o = np.array([j[0] for j in o])

    assert (o.shape == w.shape)
    return o, w


def loadSignal(_year, _region, _signal, _binning, _obs):

    observable = data_loader.sig[_year][_signal]["observables"][_obs]
    weights = data_loader.sig[_year][_signal]["weights"]
    region = data_loader.sig[_year][_signal]["regions"]

    # Select only events in a given CR
    observable, weights = filter_events(observable, weights, region, _region, _obs)

    sign, bins = np.histogram(observable, weights=weights, bins=_binning)
    data = np.maximum(sign, 0)
    error = np.sqrt(sign)

    return sign, error


def loadMC(_year, _region, _binning, _obs):
    histograms = {}
    errors = {}
 
    for bkg in BACKGROUNDS:
        
        observable = data_loader.bkg[_year][bkg]["observables"][_obs]
        weights = data_loader.bkg[_year][bkg]["weights"]
        region = data_loader.bkg[_year][bkg]["regions"]

        observable, weights = filter_events(observable, weights, region, _region, _obs)
        n, bins = np.histogram(observable, weights=weights, bins=_binning)
        histograms[bkg] = np.maximum(n, 0)

        n_err = np.sqrt(np.histogram(observable, weights=weights**2, bins=_binning)[0])
        errors[bkg] = n_err

    return histograms, errors


def regionName(_region):
    if _region == 0:
        return "Signal Region"
    else:
        return "Undefined"


def get_obs_labels(observable):
    """ Maps the ntuple key of the observable to 
    a tuple of the plot legend label and a shorthand
    for the output filename of the plot. """

    x = {
        "jetsAk4CHS/jetsAk4CHS.m_pt": (r"$p_{T}(j_{1})$ [GeV]", "Jet1PT"),
        "jetsAk4CHS/jetsAk4CHS.m_phi": (r"$\phi(j_{1})$ [rad]", "Jet1Phi"),
        "jetsAk4CHS/jetsAk4CHS.m_eta": (r"$\eta(j_{1})$", "Jet1Eta"),
        "A_mt": (r"$m_T$ of A [GeV]","AMT"),
        "H_mt": (r"$m_T$ of H [GeV]","HMT"),
        "MET": (r"missing $E_T$ [GeV]","MET"),
        "event_weight": (r"Event Weight","Weight"),
        "W_m": (r"best $W_m$ [GeV]","Wmass"),
        "HT": (r"$H_T$ [GeV]","HT"),
        "tight_b": (r"tight b tags", "tightb"), 
        }
    return x[observable.replace("Puppi", "CHS")]


def plot_bkg_vs_sig(_year, _region, _signal, _binning=np.linspace(50, 750, 70), _obs='jetsAk4CHS/jetsAk4CHS.m_pt'):

    histograms, errors = loadMC(_year, _region, _binning, _obs)
    np_histograms = np.array(list(histograms.values())) 
    
    sig, error = loadSignal(_year, _region, _signal, _binning, _obs)
    sig = sig / sum(sig) * sum(sum(np_histograms))
    error = error / sum(sig) * sum(sum(np_histograms))

    fig, axes = plt.subplots(
        nrows=2, 
        ncols=1,
        figsize=(12,10),
        gridspec_kw={"height_ratios": (4, 1)},
        sharex=True
    )

    hep.cms.label(ax=axes[0],llabel='Work in progress',data=True, lumi=YEAR_LUMI_MAP[_year], year=UL_YEAR_MAP[_year])

    np_errors = np.array(list(errors.values()))
    labels = list(histograms.keys())
    hep.histplot(np_histograms, histtype='fill', w2=np_errors, bins=_binning, stack=True, label=labels, ax=axes[0])
    hep.histplot(sig, yerr=error, histtype='errorbar', bins=_binning, markersize=13, color='k', label='Signal', ax=axes[0])

    errors_bkg = np.sqrt(sum(np_errors)**2)
    errors_sig = np.sqrt(sig)
    bkg = sum(np_histograms)
    ratio = np.nan_to_num(sig / bkg)
    errors = np.nan_to_num(ratio * np.sqrt((errors_bkg / bkg) ** 2 + (errors_sig / sig) ** 2))

    axes[0].set_yscale('log')
    axes[0].legend(ncol=2, title=f"{regionName(_region)}, {_year}", fontsize=18, 
        title_fontsize=18, frameon=True)
    axes[0].set_ylabel("Events")

    #hep.histplot(ratio, yerr=errors, histtype='errorbar', bins=_binning, stack=False, color='k', ax=axes[1])
    #axes[1].axhline(1, color='grey', linestyle='--') 
    #axes[1].set_ylabel('DATA/MC')
    #axes[1].set_ylim([0.5, 1.5])
    #axes[1].set_xlabel(get_obs_labels(_obs)[0])

    # Save Plot
    folder = os.path.join(ANALYSIS,"plots",_signal)
    if not os.path.exists(folder): os.makedirs(folder)

    image = regionName(_region).replace(' ', '') + '_' + get_obs_labels(_obs)[1] + '_' + _year

    plt.savefig(os.path.join(folder,image))
    plt.close()



def get_n_plots(obs):
    i = 0
    for year, region, signal in product(YEARS, REGIONS, SIGNALS): i += 1
    return i


if __name__ == "__main__":

    observables = {
        "jetsAk4CHS/jetsAk4CHS.m_pt": np.linspace(80, 800, 21),
        "jetsAk4CHS/jetsAk4CHS.m_phi": np.linspace(-pi, pi, 21),
        "jetsAk4CHS/jetsAk4CHS.m_eta": np.linspace(-3, 3, 21), 
        "A_mt": np.linspace(300, 2000, 51),
        "H_mt": np.linspace(300, 1500, 41),
        "MET": np.linspace(0, 800, 21),
        "event_weight": np.linspace(0,100,11),
        "W_m": np.linspace(0,400,21),
        "HT": np.linspace(100,2000,51),
        "tight_b": np.linspace(-0.25,2.25,6),
        #"electrons_tight/electrons_tight.m_pt": np.concatenate([np.linspace(0,300,20),np.linspace(350, 800, 11)]),
        #"muons_tight/muons_tight.m_pt": np.concatenate([np.linspace(0,300,20),np.linspace(350, 800, 11)]),
        #"z_pt_reco": np.linspace(0, 650, 16),
        #"a_minus_h_mass_dnn_fixedb": np.concatenate([np.linspace(80,600,20), np.linspace(700, 1300, 7)]),
        #"a_minus_h_mass_chi_sq_fixedb": np.concatenate([np.linspace(80,600,20), np.linspace(700, 1300, 7)]),
        #"h_mass_reco_dnn_fixedb": np.linspace(150, 3000, 50),
        #"a_mass_dnn_fixedb": np.linspace(250, 3500, 50),
        #"z_mass_reco": np.linspace(30, 150, 50),
        #"jetsAk4CHS/jetsAk4CHS.m_btag_DeepFlavour_probb": np.linspace(0, 1, 50),
        #"chi_sq_lowest": np.logspace(np.log(0.1), np.log(2500), 30, base=np.e),
        #"dnn_score_highest": np.linspace(0, 1, 50),
        }

    data_loader = DataLoader(observables)
    for obs, bins in observables.items():
        bar = IncrementalBar(f"Running Plots for {obs}", max=get_n_plots(obs))
        for year, region, signal in product(YEARS, REGIONS, SIGNALS):
            plot_bkg_vs_sig(year, region, signal, bins, obs)
            bar.next()
        bar.finish()

