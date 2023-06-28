#!/bin/bash

if [ "$1" == "preselection" ]; then
  echo "merge preselection"
  cd $CMSSW_BASE/src/UHH2/AZH/data/output_01_preselection/MC/UL17/
elif [ "$1" == "reconstruction" ]; then
  echo "merge reconstruction"
  cd $CMSSW_BASE/src/UHH2/AZH/data/output_02_reconstruction/MC/UL17/
else 
  echo "wrong argument"
  exit 0
fi

cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-600_MH-400_JESnominal_JERnominal_UL17.root MC.AZH_600_400_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-700_MH-450_JESnominal_JERnominal_UL17.root MC.AZH_700_450_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-750_MH-400_JESnominal_JERnominal_UL17.root MC.AZH_750_400_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-750_MH-650_JESnominal_JERnominal_UL17.root MC.AZH_750_650_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-800_MH-400_JESnominal_JERnominal_UL17.root MC.AZH_800_400_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-1000_MH-400_JESnominal_JERnominal_UL17.root MC.AZH_1000_400_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-1000_MH-850_JESnominal_JERnominal_UL17.root MC.AZH_1000_850_UL17.root

hadd MC.QCD_UL17.root *MC.QCD*.root
hadd MC.SingleTop_UL17.root *MC.ST*.root
hadd MC.TTW_UL17.root *MC.TTW*.root
hadd MC.TTZ_UL17.root *MC.TTZ*.root
hadd MC.VV_UL17.root *MC.WW* *MC.WZ* *MC.ZZ*
hadd MC.WJets_UL17.root *MC.WJets*.root
hadd MC.ZJets_UL17.root *MC.ZJets*.root *MC.DYJets*.root
hadd MC.TT_UL17.root *MC.TTTo*.root

#hadd MC.TT_UL17_hdampUp.root *TTTo*hdampUP*.root
#hadd MC.TT_UL17_hdampDown.root *TTTo*hdampDOWN*.root


cd $CMSSW_BASE/src/UHH2/AZH
