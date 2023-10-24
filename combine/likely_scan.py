#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

import os
import argparse
from fitter import execute

def log_likely_scan(workspace, outfile, r_in=0):
    r_max = max(2*r_in,2)

    HASH = workspace.replace(".root","")
    
    cmd1 = (
        f"combine -M MultiDimFit {workspace} -n .{HASH}.nominal "
        f"-m 125 --rMin -1 --rMax {r_max} --algo grid --points 30 -t -1 --setParameters r={r_in}"
    )

    cmd2 = (
        f"combine -M MultiDimFit {workspace} -n .{HASH}.snapshot "
        f"-m 125 --rMin -1 --rMax {r_max} --saveWorkspace -t -1 --setParameters r={r_in}"
    )

    cmd3 = (
        f"combine -M MultiDimFit higgsCombine.{HASH}.snapshot.MultiDimFit.mH125.root "
        f"-n .{HASH}.freezeall -m 125 --rMin -1 --rMax {r_max} --algo grid --points 30 "
        f"--freezeParameters allConstrainedNuisances --snapshotName MultiDimFit -t -1 --setParameters r={r_in}"
    )

    cmd4 = (
        f"python $CMSSW_BASE/src/CombineHarvester/CombineTools/scripts/plot1DScan.py "
        f"higgsCombine.{HASH}.nominal.MultiDimFit.mH125.root --others "
        f"'higgsCombine.{HASH}.freezeall.MultiDimFit.mH125.root:FreezeAll:{r_max}' -o {outfile}"
    )

    cmd5 = (
        f"rm higgsCombine.{HASH}*"
    )

    execute(cmd1,cmd2,cmd3,cmd4,cmd5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("masses")
    args = parser.parse_args()
    masses = args.masses.replace("MA-","").replace("MH-","")

    workspace = f"WS_{masses}_MET_inv_all.root"
    outfile = workspace.replace("WS","SCAN").replace(".root","")
    log_likely_scan(workspace,outfile)

#combine -M MultiDimFit -n .nominal --algo grid --points 50 --redefineSignalPOIs r -P CMS_toppt_b --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 -m 125 -d workspace.root ..

