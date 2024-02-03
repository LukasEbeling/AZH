#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import argparse
from plot_utils import REGIONS
from combine import Combine
import sys

sys.path.append("../factory")

TEMPLATES = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/data/output_03_templates"


if __name__ == "__main__":
    parser = argparse.ArgumentParser() #eg ./run_impatcs.py 1000_400 MET 0.035
    parser.add_argument("masses")
    parser.add_argument("var")
    parser.add_argument("r")
    args = parser.parse_args()
    masses = args.masses; var = args.var; r = args.r   
    
    masses = masses.replace("MA-","").replace("MH-","")

    cards = [TEMPLATES+f"/UL17/AZH_{masses}_{var}_inv_{region}.dat" for region in REGIONS.keys()]
    outcard = TEMPLATES+f"/UL17/AZH_{masses}_{var}_inv_all.dat"
    workspace = f"{masses}_{var}_inv_all.root"
    
    combine = Combine(initilize=True)
    combine.combine_cards(cards,outcard)
    combine.create_workspace(outcard,workspace)
    combine.calc_impacts(workspace,r)