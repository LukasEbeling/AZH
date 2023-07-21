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

hadd -f MC.QCD_UL17.root $(ls -S *MC.QCD*) 
hadd -f MC.SingleTop_UL17.root $(ls -S *MC.ST*)
hadd -f MC.TTW_UL17.root $(ls -S *MC.TTW*)
hadd -f MC.TTZ_UL17.root $(ls -S *MC.TTZ*)
hadd -f MC.VV_UL17.root $(ls -S *MC.WW*) $(ls -S *MC.WZ*) $(ls -S *MC.ZZ*) 
hadd -f MC.WJets_UL17.root $(ls -S *MC.WJets*)
hadd -f MC.ZJets_UL17.root $(ls -S *MC.ZJets*) $(ls -S *MC.DYJets*)
hadd -f MC.TT_UL17.root $(ls -S *MC.TTTo*)



cd $CMSSW_BASE/src/UHH2/AZH
