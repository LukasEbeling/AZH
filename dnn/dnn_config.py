#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import os

#CMSSW_BASE = os.environ.get('CMSSW_BASE')
CMSSW_BASE = '/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/'
CACHE = os.path.join(CMSSW_BASE, "src/UHH2/AZH/data/output_03_templates/cache")

### dnn structure ###
NODES = 512
EPOCHS = 10
BATCHSIZE = 320

### targets for bkg processes ###
TARGETS = {
    #Top node
    "TTZ": [0,1,0,0,0],
    "TTW": [0,1,0,0,0],
    "SingleTop": [0,1,0,0,0],
    "TT": [0,1,0,0,0],

    #WJets + VV node
    "WJets_ljet": [0,0,1,0,0],
    "WJets_bjet": [0,0,1,0,0], 
    "VV": [0,0,1,0,0],

    #DYJets node
    "DYJets_ljet": [0,0,0,1,0],  
    "DYJets_bjet": [0,0,0,1,0],
    
    #QCD node
    "QCD": [0,0,0,0,1],
}

### parquet files to load ###
LOAD = [
    "jetsAk4CHS.m_pt_1",
    "jetsAk4CHS.m_pt_2",
    "jetsAk4CHS.m_pt_3",
    "jetsAk4CHS.m_eta_1",
    "jetsAk4CHS.m_eta_2",
    "jetsAk4CHS.m_eta_3",
    "jetsAk4CHS.m_phi_1",
    "jetsAk4CHS.m_phi_2",
    "jetsAk4CHS.m_phi_3",
    "jetsAk4CHS.m_btag_DeepFlavour_probbb_1",
    "jetsAk4CHS.m_btag_DeepFlavour_probbb_2",
    "jetsAk4CHS.m_btag_DeepFlavour_probbb_3",
    "jetsAk4CHS.m_btag_DeepFlavour_probb_1",
    "jetsAk4CHS.m_btag_DeepFlavour_probb_2",
    "jetsAk4CHS.m_btag_DeepFlavour_probb_3",
    "jetsAk4CHS.m_btag_DeepFlavour_problepb_1",
    "jetsAk4CHS.m_btag_DeepFlavour_problepb_2",
    "jetsAk4CHS.m_btag_DeepFlavour_problepb_3",
    "MET",
    "m_H",
    "mt_A",
    "delta_phi",
    "delta_eta",
    "num_b",
    "num_l",
    "num_j",
    "HT",
    "event_weight",
    "kfold",
    "region",
]

### input features - observables considered by dnn ###
OBSERVABLES = [
    "jetsAk4CHS.m_pt_1",
    "jetsAk4CHS.m_pt_2",
    "jetsAk4CHS.m_pt_3",
    "jetsAk4CHS.m_eta_1",
    "jetsAk4CHS.m_eta_2",
    "jetsAk4CHS.m_eta_3",
    "MET",
    "m_H",
    "mt_A",
    "delta_phi",
    "delta_eta",
    "num_b",
    "num_l",
    "num_j",
    "HT",
    "bscore1",
    "bscore2",
    "bscore3",
    "delta_phi_12",
    "delta_phi_13",
    "delta_phi_23",
]

