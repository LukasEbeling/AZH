#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import argparse
from plot_utils import REGIONS
from combine import Combine
import sys

sys.path.append("../factory")

TEMPLATES = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/data/output_03_templates"# combination with lepton channel

if __name__ == "__main__":
    parser = argparse.ArgumentParser() #eg ./run_combination.py 1000_400 MET
    parser.add_argument("masses")
    parser.add_argument("var")
    args = parser.parse_args()
    masses = args.masses; var = args.var   
    
    masses = masses.replace("MA-","").replace("MH-","")
    
    lep_regions = [
        'SR-6j-2b',
        'CR-6j-2b-SB',
        'CR-5j-0b',
        'SR-5j-1b',
        'CR-6j-0b',
        'CR-6j-1b-SB',
        'SR-5j-2b',
        'SR-6j-1b',
        'CR-5j-2b-SB',
        'CR-5j-1b-SB',
    ]

    elec = [TEMPLATES+f"/UL17/AZH_{masses}_2DEllipses_diElectron_{region}.dat" for region in lep_regions]
    muon = [TEMPLATES+f"/UL17/AZH_{masses}_2DEllipses_diMuon_{region}.dat" for region in lep_regions]
    outcard = TEMPLATES+f"/UL17/AZH_{masses}_2DEllipses_lep_all.dat"
    workspace = f"{masses}_2DEllipses_lep_all.root"
    logfile = f"AZH_{masses}_2DEllipses_lep_all.log"

    combine = Combine(initilize=True)
    combine.combine_cards(elec+muon,outcard)
    combine.create_workspace(outcard,workspace)
    combine.run_blind(workspace,logfile)
    
    inv = [TEMPLATES+f"/UL17/AZH_{masses}_{var}_inv_{region}.dat" for region in REGIONS.keys()]
    outcard = TEMPLATES+f"/UL17/AZH_{masses}_{var}_combined.dat"
    workspace = f"{masses}_{var}_combined_all.root"
    logfile = f"AZH_{masses}_{var}_combined_all.log"
    combine.combine_cards(elec+muon+inv,outcard)
    combine.create_workspace(outcard,workspace)
    combine.run_blind(workspace,logfile)
    
