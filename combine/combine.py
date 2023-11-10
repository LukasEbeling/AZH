import os

# execute some commands after initializing combine
def execute(*cmds):
    init_old = (
        "source /cvmfs/cms.cern.ch/cmsset_default.sh; "
        "source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh; "
        "cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit; "
        "cmsenv; "
        "cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/combine/; "
    )
    init = "source /afs/desy.de/user/e/ebelingl/combine_ini; "
    commands = [cmd+"; " for cmd in cmds]
    joined = "".join(commands)
    os.system(init+joined)


### interface with combine tool ###
class Combine():

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

    # run fit diagnostics
    def fit(card,workspace,output): 
        short = card.replace(".dat","").replace("/","_")
        execute(
            f"combine -M FitDiagnostics -d {workspace} "
            f"-n _{short} -m 125 -v 1 --rMin -1 --rMax 1 "
            "--expectSignal 0.0 -t -1 --saveWithUncertainties --cminDefaultMinimizerStrategy 0 --saveShapes "
        )
        execute(
            f"PostFitShapesFromWorkspace -d {card} -w {workspace} "
            f"--output {output} -m 125 "
            "--sampling --print --covariance --total-shapes "
            f"--postfit -f fitDiagnostics_{short}.root:fit_s"
        )
        execute(f"rm fitDiagnostics_{short}.root")
        execute(f"rm higgsCombine_{short}*")


    # perform a log likelyhood scan of poi
    def log_likely_scan(workspace, outfile, r_in=0):
        r_max = max(2*r_in,2)

        HASH = workspace.replace(".root","")
        
        execute(
            f"combine -M MultiDimFit {workspace} -n .{HASH}.nominal "
            f"-m 125 --rMin -1 --rMax {r_max} --algo grid --points 30 -t -1 --setParameters r={r_in}"
        )

        execute(
            f"combine -M MultiDimFit {workspace} -n .{HASH}.snapshot "
            f"-m 125 --rMin -1 --rMax {r_max} --saveWorkspace -t -1 --setParameters r={r_in}"
        )

        execute(
            f"combine -M MultiDimFit higgsCombine.{HASH}.snapshot.MultiDimFit.mH125.root "
            f"-n .{HASH}.freezeall -m 125 --rMin -1 --rMax {r_max} --algo grid --points 30 "
            f"--freezeParameters allConstrainedNuisances --snapshotName MultiDimFit -t -1 --setParameters r={r_in}"
        )

        execute(
            f"python $CMSSW_BASE/src/CombineHarvester/CombineTools/scripts/plot1DScan.py "
            f"higgsCombine.{HASH}.nominal.MultiDimFit.mH125.root --others "
            f"'higgsCombine.{HASH}.freezeall.MultiDimFit.mH125.root:FreezeAll:{r_max}' -o {outfile}"
        )

        execute(
            f"rm higgsCombine.{HASH}*"
        )