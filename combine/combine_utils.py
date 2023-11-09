import os

### interface with combine tool ###
### execute all methods in shell where combine was initialized ###
class Combine():
    # combine a list of given datacards
    def combine_cards(card_list,card_name):
        cmd = f"combineCards.py "
        for i, card in enumerate(card_list):
            cmd += f"Name{i}={card} "
        cmd += f"> {card_name} "
        os.system(cmd)

    # create workspace from datacard
    def create_workspace(card,workspace):
        if os.path.exists(workspace): return
        cmd = f"text2workspace.py {card} -m 125 -o {workspace}"
        os.system(cmd)

    # calculate expected limit from workspace
    def run_blind(workspace, logfile, option=""):
        cmd = f"combine -M AsymptoticLimits -m 125 --run blind -d {workspace} {option} --rMin -1 --rMax 1| grep Expected &> {logfile}"
        os.system(cmd)

    # run fit diagnostics
    def fit(card,workspace,output): 
        short = card.replace(".dat","").replace("/","_")
        os.system(
            f"combine -M FitDiagnostics -d {workspace} "
            f"-n _{short} -m 125 -v 1 --rMin -1 --rMax 1 "
            "--expectSignal 0.0 -t -1 --saveWithUncertainties --cminDefaultMinimizerStrategy 0 --saveShapes "
        )
        os.system(
            f"PostFitShapesFromWorkspace -d {card} -w {workspace} "
            f"--output {output} -m 125 "
            "--sampling --print --covariance --total-shapes "
            f"--postfit -f fitDiagnostics_{short}.root:fit_s"
        )
        os.system(f"rm fitDiagnostics_{short}.root")
        os.system(f"rm higgsCombine_{short}*")


    # perform a log likelyhood scan of poi
    def log_likely_scan(workspace, outfile, r_in=0):
    r_max = max(2*r_in,2)

    HASH = workspace.replace(".root","")
    
    os.system(
        f"combine -M MultiDimFit {workspace} -n .{HASH}.nominal "
        f"-m 125 --rMin -1 --rMax {r_max} --algo grid --points 30 -t -1 --setParameters r={r_in}"
    )

    os.system(
        f"combine -M MultiDimFit {workspace} -n .{HASH}.snapshot "
        f"-m 125 --rMin -1 --rMax {r_max} --saveWorkspace -t -1 --setParameters r={r_in}"
    )

    os.system(
        f"combine -M MultiDimFit higgsCombine.{HASH}.snapshot.MultiDimFit.mH125.root "
        f"-n .{HASH}.freezeall -m 125 --rMin -1 --rMax {r_max} --algo grid --points 30 "
        f"--freezeParameters allConstrainedNuisances --snapshotName MultiDimFit -t -1 --setParameters r={r_in}"
    )

   os.system(
        f"python $CMSSW_BASE/src/CombineHarvester/CombineTools/scripts/plot1DScan.py "
        f"higgsCombine.{HASH}.nominal.MultiDimFit.mH125.root --others "
        f"'higgsCombine.{HASH}.freezeall.MultiDimFit.mH125.root:FreezeAll:{r_max}' -o {outfile}"
    )

    os.system(
        f"rm higgsCombine.{HASH}*"
    )