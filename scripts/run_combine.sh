cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit
cmsenv
. env_standalone.sh
make -j 4
cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/combine/UL17/
    
    
combine_cards() {
    echo combining cards for $1 
    combineCards.py Name1=INV_${1}_met_inv_SignalRegion.dat Name2=INV_${1}_met_inv_CR_0B.dat Name3=INV_${1}_met_inv_CR_1L.dat Name4=INV_${1}_met_inv_CR_lowmet.dat > INV_${1}_met.dat
    #combineCards.py Name1=INV_${1}_Amt_inv_SignalRegion.dat Name2=INV_${1}_Amt_inv_CR_0B.dat > INV_${1}_Amt.dat
}

calculate_limit() {
    echo AsymptoticLimits $1 processing...
    combine -M AsymptoticLimits -m 125 --run blind -d ${path}INV_${1}_met_inv_SignalRegion.dat --rMin -1 --rMax 1| grep Expected &> ${path}AsymptoticLimits_AZHToInv_${1}_met_SR.log
    combine -M AsymptoticLimits -m 125 --run blind -d ${path}INV_${1}_met.dat --rMin -1 --rMax 1| grep Expected &> ${path}AsymptoticLimits_AZHToInv_${1}_met.log
    #combine -M AsymptoticLimits -m 125 --run blind -d ${path}INV_${1}_Amt.dat --rMin -1 --rMax 1| grep Expected &> ${path}AsymptoticLimits_AZHToInv_${1}_Amt.log

}

create_workspaces() {
    echo text2workspace.py $1
    text2workspace.py INV_${1}_met.dat -m 125 -o INV_${1}_met_workspace.root
    #text2workspace.py INV_${1}_Amt.dat -m 125 -o INV_${1}_Amt_workspace.root
}

plot_impacts() {
    echo plotting impacts for $1
    combineTool.py -M Impacts -d INV_${1}_met_workspace.root -m 125 --doInitialFit --robustFit 1 --rMin -1 --rMax 1 --expectSignal 0.05 -t -1
    combineTool.py -M Impacts -d INV_${1}_met_workspace.root -m 125 --doFits --robustFit 1 --rMin -1 --rMax 1 --expectSignal 0.05 -t -1
    combineTool.py -M Impacts -d INV_${1}_met_workspace.root -m 125 -o impacts.json --rMin -1 --rMax 1 --expectSignal 0.05 -t -1
    plotImpacts.py -i impacts.json -o impacts_${1}_met

    #combineTool.py -M Impacts -d INV_${1}_Amt_workspace.root -m 125 --doInitialFit --robustFit 1 --rMin -1 --rMax 1 --expectSignal 0.05 -t -1
    #combineTool.py -M Impacts -d INV_${1}_Amt_workspace.root -m 125 --doFits --robustFit 1 --rMin -1 --rMax 1 --expectSignal 0.05 -t -1
    #combineTool.py -M Impacts -d INV_${1}_Amt_workspace.root -m 125 -o impacts.json --rMin -1 --rMax 1 --expectSignal 0.05 -t -1
    #plotImpacts.py -i impacts.json -o impacts_${1}_Amt

}


masses=("1000_400" "1000_850" "600_400" "700_450" "750_400" "750_650" "800_400")

if [ "$1" == "limits" ]; then
    for m in ${masses[@]}; do combine_cards $m; done
    for m in ${masses[@]}; do calculate_limit $m; done
elif [ "$1" == "impacts" ]; then
    for m in ${masses[@]}; do combine_cards $m; done
    for m in ${masses[@]}; do create_workspaces $m; done
    for m in ${masses[@]}; do plot_impacts $m; done
    rm combine_logger.out
    rm impacts.json
else 
  echo "wrong argument"
  exit 0
fi

rm higgsCombine*