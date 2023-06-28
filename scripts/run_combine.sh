#!/bin/bash
path="/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/combine/UL17/"

preparation() {
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
    cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit
    cmsenv
    . env_standalone.sh
    make -j 4
    cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/combine/UL17/
}


combine_cards() {
    local mass=$1
    local var=$2

    echo combining cards for $mass
    combineCards.py Name1=AZH_${mass}_${var}_inv_SignalRegion.dat \
        Name2=AZH_${mass}_${var}_inv_CR_0B.dat Name3=AZH_${mass}_met_inv_CR_1L.dat \
        Name4=AZH_${mass}_${var}_inv_CR_lowmet.dat > AZH_${mass}_${var}.dat
}


calculate_limit() {
    local mass=$1
    local var=$2

    echo AsymptoticLimits for $mass
    combine -M AsymptoticLimits -m 125 --run blind -d ${path}AZH_${mass}_${var}_inv_SignalRegion.dat \
        --rMin -1 --rMax 1| grep Expected &> ${path}expected_${mass}_${var}_SR.log
    combine -M AsymptoticLimits -m 125 --run blind -d ${path}AZH_${mass}_${var}.dat \
        --rMin -1 --rMax 1| grep Expected &> ${path}expected_${mass}_${var}.log
}


create_workspaces() {
    local mass=$1
    local var=$2

    echo text2workspace.py for ${mass}
    text2workspace.py AZH_${mass}_${var}.dat -m 125 -o AZH_${mass}_${var}_workspace.root
}


plot_impacts() {
    local mass=$1
    local var=$2
    
    echo plotting impacts for ${mass}
    combineTool.py -M Impacts -d AZH_${mass}_${var}_workspace.root -m 125 \
        --doInitialFit --robustFit 1 --rMin -1 --rMax 1 --expectSignal 0.05 -t -1
    combineTool.py -M Impacts -d AZH_${mass}_${var}_workspace.root -m 125 \
        --doFits --robustFit 1 --rMin -1 --rMax 1 --expectSignal 0.05 -t -1
    combineTool.py -M Impacts -d AZH_${mass}_${var}_workspace.root -m 125 -o impacts_${mass}_${var}.json \
        --rMin -1 --rMax 1 --expectSignal 0.05 -t -1
    plotImpacts.py -i impacts_${mass}_${var}.json -o impacts_${mass}_${var}
}



preparation
combine_cards $1 $2
calculate_limit $1 $2

if [ "$3" = "impacts" ]; then
    create_workspaces $1 $2
    plot_impacts $1 $2
fi
