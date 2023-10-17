#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
from collections import OrderedDict as od
from itertools import product
from math import pi, sqrt
import os 

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pandas as pd
from progress.bar import IncrementalBar
import uproot


plt.style.use(hep.style.CMS)

#Plotting the number of leptons @gen and @reco 
#Relevant for TT -> semileptonic

#Config
ANALYSIS = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/"
PRESELECTION = ANALYSIS + "data/output_01_preselection/MC/UL17/"
FILE = PRESELECTION + "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_JESnominal_JERnominal_UL17.root"


OBSERVABLES = {
    "lep": "leptons",
    "id": "GenParticles.m_pdgId",
    "pt": "GenParticles.m_pt",
    "eta": "GenParticles.m_eta",
    "phi": "GenParticles.m_phi",
    "j_pt": "jetsAk4CHS.m_pt",
    "j_eta": "jetsAk4CHS.m_eta",
    "j_phi": "jetsAk4CHS.m_phi",
}

class object:
    pt: float = 0
    eta: float = 0
    phi: float = 0

    def __init__(self,pt,eta,phi):
        self.pt = pt
        self.eta = eta 
        self.phi = phi

class jet(object): pass
class electron(object): pass
class muon(object): pass

class event:
    electrons: list = []
    muons: list = []
    jets: list = []

    def __init__(self,ids,gen_pts,gen_etas,gen_phis,jet_pts,jet_etas,jet_phis):
        self.electrons = []
        self.muons = []
        self.jets = []

        self.create_jets(jet_pts,jet_etas,jet_phis)
        self.create_gen(ids,gen_pts,gen_etas,gen_phis)
        self.clean_collection(self.electrons)
        self.clean_collection(self.muons)
        print("new event")

    def create_jets(self,jet_pts,jet_etas,jet_phis):
        self.jets = []
        for i in range(len(jet_pts)):
            if jet_pts[i] < 30: continue
            #if abs(jet_etas[i]) > 2.4: continue
            new_jet = jet(jet_pts[i],jet_etas[i],jet_phis[i])
            self.jets.append(new_jet)

    def create_gen(self,ids,gen_pts,gen_etas,gen_phis):
        for i in range(len(ids)):
            
            if gen_pts[i] < 20: continue
            if gen_etas[i] > 2.4: continue
            
            if abs(ids[i]) == 11:
                new_elec = electron(gen_pts[i],gen_etas[i],gen_phis[i])
                self.electrons.append(new_elec)
            elif abs(ids[i]) == 13:
                new_muon = muon(gen_pts[i],gen_etas[i],gen_phis[i])
                self.muons.append(new_muon)
            else: continue

    def clean_collection(self,particles):
        if len(self.jets) == 0: return True
        for p in particles:
            if self.min_distance(p) < 0.8: particles.remove(p)

    def min_distance(self,p):
        distances = [self.distance(p,j) for j in self.jets]
        return min(distances)

    def distance(self,a,b):
        delta_phi = abs(a.phi-b.phi)
        if delta_phi > pi: delta_phi = 2*pi - delta_phi
        delta_eta = abs(a.eta-b.eta)
        return sqrt(delta_phi**2+delta_eta**2)

    def get_gen(self):
        return len(self.electrons)+len(self.muons)




class dataLoader():

    data: dir = {}
    events: list = []

    def __init__(self):
        self.load_observables()
        self.build_events()
    
    def load_observables(self):
        for key in OBSERVABLES.keys():
            branch = f"AnalysisTree/{OBSERVABLES[key]}" 
            self.data[key] = self.load(branch)      
        
    def load(self, branch): 
        with uproot.open(FILE) as file:
            return file[branch].array(library="np")
        
    def build_events(self):
        g_ids = self.data["id"]
        g_pts = self.data["pt"]
        g_etas = self.data["eta"]
        g_phis = self.data["phi"]
        j_pts = self.data["j_pt"]
        j_etas = self.data["j_eta"]
        j_phis = self.data["j_phi"]

        for i in range(len(g_ids)):
            new_event = event(g_ids[i],g_pts[i],g_etas[i],g_phis[i],j_pts[i],j_etas[i],j_phis[i])
            self.events.append(new_event)

    def get_gen(self):
        gen_leptons = [event.get_gen() for event in self.events]
        #return sum(gen_leptons)/len(gen_leptons)
        return gen_leptons
    
    def get_reco(self):
        #return sum(self.data["lep"])/len(self.data["lep"])
        return self.data["lep"]
    
def plot(num_gen,num_reco):        
    plt.hist(num_gen,bins=[0,1,2,3,4],histtype="step",color="red",label="gen")
    plt.hist(num_reco,bins=[0,1,2,3,4],histtype="step",label="reco")
    plt.xticks([0.5,1.5,2.5,3.5,4.5],['0', '1', '2','3','4'])
    plt.legend(loc='upper right', fontsize='small')
    plt.savefig(os.path.join(ANALYSIS,"genvsreco"))



if __name__ == "__main__":
    data_loader = dataLoader()
    gen = data_loader.get_gen()
    reco = data_loader.get_reco()

    plot(gen,reco)