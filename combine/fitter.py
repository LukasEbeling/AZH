#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import os
import argparse
from itertools import product

INIT = "/afs/desy.de/user/e/ebelingl/combine_ini"
ANALYSIS = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/combine"
YEAR = "UL17" # later via argparse


INV_REGIONS = [
    "SR_6J",
    "SR_5J",
    "IR_1B_5J",
    "IR_1B_6J",
    "IR_0B_5J",
    "IR_0B_6J",
    "LR_2B_5J",
    "LR_2B_6J",
    "LR_1B_5J",
    "LR_1B_6J",
    "LR_0B_5J",
    "LR_0B_6J",
]

INV_CHANNELS = [
    "inv"
]

LEP_REGIONS = [
    "CR-5j-0b",
    "CR-5j-1b-SB",
    "CR-5j-2b-SB",
    "CR-6j-0b",
    "CR-6j-1b-SB",
    "CR-6j-2b-SB",
    "SR-5j-1b",
    "SR-5j-2b",
    "SR-6j-1b",
    "SR-6j-2b",
]

LEP_CHANNELS = [
    "diMuon",
    "diElectron"
]

def execute(command):
    cmd = f"source {INIT}; "
    cmd += f"cd {ANALYSIS}/{YEAR}; "
    cmd += command
    os.system(cmd)

def clean_up():
    os.system(f"rm {ANALYSIS}/{YEAR}/higgsCombine*")
    os.system(f"rm {ANALYSIS}/{YEAR}/combine_logger.out")

def combine_cards(card_list,card_name):
    cmd = f"combineCards.py "
    for i, card in enumerate(card_list):
        cmd += f"Name{i}={card}.dat "
    cmd += f"> {card_name}.dat"
    execute(cmd)

def run_blind(card):
    cmd = f"combine -M AsymptoticLimits -m 125 --run blind -d {card}.dat --rMin -1 --rMax 1| grep Expected &> {card}.log"
    execute(cmd)

def create_workspace(card, workspace): 
    cmd = f"text2workspace.py {card}.dat -m 125 -o {workspace}.root"
    execute(cmd)

def plot_impacts(card, pdf):
    workspace = card + "_workspace"
    create_workspace(card, workspace)
    cmd = f"combineTool.py -M Impacts -d {workspace}.root -m 125 --doInitialFit --robustFit 1 --rMin -1 --rMax 1 --expectSignal 0.05 -t -1; "
    cmd += f"combineTool.py -M Impacts -d {workspace}.root -m 125 --doFits --robustFit 1 --rMin -1 --rMax 1 --expectSignal 0.05 -t -1; "
    cmd += f"combineTool.py -M Impacts -d {workspace}.root -m 125 -o {card}_impacts.json --rMin -1 --rMax 1 --expectSignal 0.05 -t -1; "
    cmd += f"plotImpacts.py -i {card}_impacts.json -o {pdf}"   
    execute(cmd)


def run_invisible(mA,mH,var,impacts=False):
    cards = [f"AZH_{mA}_{mH}_{var}_{ch}_{region}" for ch,region in product(INV_CHANNELS, INV_REGIONS)]
    outcard = f"INV_{mA}_{mH}_{var}_UL17"
    combine_cards(cards, outcard)
    run_blind(outcard)
    if impacts: plot_impacts(outcard, f"impacts_{mA}_{mH}_{var}")

def run_leptonic(mA,mH,var="2DEllipses"):
    cards = [f"AZH_MA-{mA}_MH-{mH}_{var}_{ch}_{reg}" for ch, reg in product(LEP_CHANNELS, LEP_REGIONS)]
    outcard = f"LEP_{mA}_{mH}_{var}_UL17"
    combine_cards(cards, outcard)
    run_blind(outcard)

def run_both(mA,mH,var1,var2="2DEllipses"):
    cards = [f"AZH_MA-{mA}_MH-{mH}_{var2}_{ch}_{reg}" for ch, reg in product(LEP_CHANNELS, LEP_REGIONS)]
    cards += [f"AZH_{mA}_{mH}_{var1}_{ch}_{region}" for ch,region in product(INV_CHANNELS, INV_REGIONS)]
    outcard = f"combined_{mA}_{mH}_{var}_UL17"
    combine_cards(cards, outcard)
    run_blind(outcard)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mA")
    parser.add_argument("mH")
    parser.add_argument("var")
    parser.add_argument('--impacts', action='store_true')
    args = parser.parse_args()
    mA = args.mA; mH = args.mH; var = args.var; impacts = args.impacts

    run_invisible(mA,mH,var,impacts)
    run_leptonic(mA,mH)
    run_both(mA,mH,var)
    clean_up()

    