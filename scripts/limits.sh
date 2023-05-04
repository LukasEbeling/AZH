path=/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/combine/UL17/

# AsymptoticLimits
run_asymptotic_limits() {
    echo AsymptoticLimits $1 processing...
    echo ${path}INV_${1}_Amt_inv_SignalRegion.dat
    combine -M AsymptoticLimits -m 125 --run blind -d ${path}INV_${1}_Amt_inv_SignalRegion.dat | grep Expected &> ${path}AsymptoticLimits_AZHToInv_${1}_Amt.log
    echo ${path}INV_${1}_Hmt_inv_SignalRegion.dat
    combine -M AsymptoticLimits -m 125 --run blind -d ${path}INV_${1}_Hmt_inv_SignalRegion.dat | grep Expected &> ${path}AsymptoticLimits_AZHToInv_${1}_Hmt.log
    echo ${path}INV_${1}_met_inv_SignalRegion.dat
    combine -M AsymptoticLimits -m 125 --run blind -d ${path}INV_${1}_met_inv_SignalRegion.dat | grep Expected &> ${path}AsymptoticLimits_AZHToInv_${1}_met.log
    echo ${path}INV_${1}_2DEllipses_inv_SignalRegion.dat
    combine -M AsymptoticLimits -m 125 --run blind -d ${path}INV_${1}_2DEllipses_inv_SignalRegion.dat | grep Expected &> ${path}AsymptoticLimits_AZHToInv_${1}_2DEllipses.log

}

cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit
cmsenv
. env_standalone.sh
make -j 4

cd $path
masses=("1000_400" "1000_850" "600_400" "700_450" "750_400" "750_650" "800_400")

for m in ${masses[@]}; do run_asymptotic_limits $m;done