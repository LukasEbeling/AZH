#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import os
from itertools import product

REGIONS = [
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

LEP_MASSES = [
    "600_400"
    "800_400"
    "750_400"
    "600_400"
]

# execute an arbitrary number of commands after sourcing combine init file
def execute(*args):
    init = "source /afs/desy.de/user/e/ebelingl/combine_ini; "
    cd = "cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/combine/UL17; "
    commands = "".join([command +"; " for command in args])
    os.system(init+cd+commands)
    #print("TEST: "+init+cd+commands+"\n")

# combine a list of given datacards
def combine_cards(card_list,card_name):
    cmd = f"combineCards.py "
    for i, card in enumerate(card_list):
        cmd += f"Name{i}={card} "
    cmd += f"> {card_name} "
    execute(cmd)

# create workspace from datacard
def create_workspace(card,workspace):
    if os.path.exists(workspace): return
    cmd = f"text2workspace.py {card} -m 125 -o {workspace}"
    execute(cmd)

# calculate expected limit from workspace
def run_blind(workspace, logfile, option=""):
    cmd = f"combine -M AsymptoticLimits -m 125 --run blind -d {workspace} {option} --rMin -1 --rMax 1| grep Expected &> {logfile}"
    execute(cmd)





    