#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import argparse
from utils import REGION_ID_MAP,Combine 


if __name__ == "__main__":
    parser = argparse.ArgumentParser() #eg ./simple_limits 1000_400 MET
    parser.add_argument("masses")
    parser.add_argument("var")
    args = parser.parse_args()
    masses = args.masses; var = args.var   
    
    masses = masses.replace("MA-","").replace("MH-","")

    cards = [f"UL17/AZH_{masses}_{var}_inv_{region}.dat" for region in REGION_ID_MAP.keys()]
    outcard = f"UL17/AZH_{masses}_{var}_inv_all.dat"
    workspace = f"UL17/WS_{masses}_{var}_inv_all.root"
    logfile = outcard.replace(".dat",".log")
    
    Combine.combine_cards(cards,outcard)
    Combine.create_workspace(outcard,workspace)
    Combine.run_blind(workspace,logfile)