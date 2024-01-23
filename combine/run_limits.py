#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import argparse
from plot_utils import REGIONS
from combine import Combine
import sys

sys.path.append("../factory")

TEMPLATES = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/data/output_03_templates"


if __name__ == "__main__":
    parser = argparse.ArgumentParser() #eg ./run_limits.py 1000_400 MET
    parser.add_argument("masses")
    parser.add_argument("var")
    args = parser.parse_args()
    masses = args.masses; var = args.var   
    
    masses = masses.replace("MA-","").replace("MH-","")

    cards = [TEMPLATES+f"/UL17/AZH_{masses}_{var}_inv_{region}.dat" for region in REGIONS.keys()]
    outcard = TEMPLATES+f"/UL17/AZH_{masses}_{var}_inv_all.dat"
    workspace = f"{masses}_{var}_inv_all.root"
    logfile = f"AZH_{masses}_{var}_inv_all.log"
    
    combine = Combine(initilize=True)
    #combine = Combine()
    combine.combine_cards(cards,outcard)
    combine.create_workspace(outcard,workspace)
    combine.run_blind(workspace,logfile)


    # individual signal regions
    for region in ['SR_2B_6J','SR_2B_5J','SR_1B_6J','SR_1B_5J']:
        card = TEMPLATES+f"/UL17/AZH_{masses}_{var}_inv_{region}.dat"
        workspace = f"AZH_{masses}_{var}_inv_{region}.root"
        logfile = f"AZH_{masses}_{var}_inv_{region}.log"
        combine.create_workspace(card,workspace)
        combine.run_blind(workspace,logfile)