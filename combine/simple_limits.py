#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import fitter
import argparse
from fitter import REGIONS,create_workspace,combine_cards,run_blind 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("masses")
    parser.add_argument("var")
    args = parser.parse_args()
    masses = args.masses; var = args.var   
    
    masses.replace("MA-","").replace("MH-","")

    cards = [f"AZH_{masses}_{var}_inv_{region}.dat" for region in REGIONS]
    outcard = f"AZH_{masses}_{var}_inv_all.dat"
    workspace = f"WS_{masses}_{var}_inv_all.root"
    logfile = outcard.replace(".dat",".log")
    
    combine_cards(cards,outcard)
    create_workspace(outcard,workspace)
    run_blind(workspace,logfile)