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
    workspace = f"tmp/{masses}_{var}_inv_all.root"
    logfile = f"tmp/AZH_{masses}_{var}_inv_all.log"
    
    Combine.combine_cards(cards,outcard)
    Combine.create_workspace(outcard,workspace)
    Combine.run_blind(workspace,logfile)