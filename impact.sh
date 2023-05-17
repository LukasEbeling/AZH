cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit
cmsenv
. env_standalone.sh
make -j 4

cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/combine/UL17

text2workspace.py INV_1000_850_Amt_inv_SignalRegion.dat -m 125

combineTool.py -M Impacts -d INV_1000_850_Amt_inv_SignalRegion.dat.root -m 125 --doInitialFit --robustFit 1 --rMin -1 --rMax 1 --expectSignal 0.5 -t -1

combineTool.py -M Impacts -d INV_1000_850_Amt_inv_SignalRegion.dat.root -m 125 --robustFit 1 --doFits --rMin -1 --rMax 1 --expectSignal 0.5 --parallel 10 -t -1

combineTool.py -M Impacts -d INV_1000_850_Amt_inv_SignalRegion.dat.root -m 125 -o impacts.json --rMin -1 --rMax 1 --expectSignal 0.5 -t -1

plotImpacts.py -i impacts.json -o impacts

rm INV_1000_850_Amt_inv_SignalRegion.dat.root
rm higgsCombine*
rm combine_logger.out
rm impacts.json

cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH