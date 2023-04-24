#!/bin/bash

if [ "$1" == "preselection" ]; then
  echo "merge preselection"
  cd $CMSSW_BASE/src/UHH2/AZH/data/output_01_preselection/MC/UL17/
elif [ "$1" == "reconstruction" ]; then
  echo "merge reconstruction"
  cd $CMSSW_BASE/src/UHH2/AZH/data/output_02_reconstruction/MC/UL17/
else 
  echo "wrong argument"
fi

cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-600_MH-400_JESnominal_JERnominal_UL17.root MC.INV_600_400_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-700_MH-450_JESnominal_JERnominal_UL17.root MC.INV_700_450_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-750_MH-400_JESnominal_JERnominal_UL17.root MC.INV_750_400_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-750_MH-650_JESnominal_JERnominal_UL17.root MC.INV_750_650_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-800_MH-400_JESnominal_JERnominal_UL17.root MC.INV_800_400_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-1000_MH-400_JESnominal_JERnominal_UL17.root MC.INV_1000_400_UL17.root
cp uhh2.AnalysisModuleRunner.MC.AToZHToInv_MA-1000_MH-850_JESnominal_JERnominal_UL17.root MC.INV_1000_850_UL17.root

cp uhh2.AnalysisModuleRunner.MC.QCD_HT300to500_JESnominal_JERnominal_UL17.root MC.HT300to500_UL17.root
cp uhh2.AnalysisModuleRunner.MC.QCD_HT500to700_JESnominal_JERnominal_UL17.root MC.HT500to700_UL17.root

hadd MC.HT700toInf_UL17.root *MC.QCD_HT700* *MC.QCD_HT1000* *MC.QCD_HT1500* *MC.QCD_HT2000* 

hadd MC.QCD_UL17.root *MC.QCD*.root
hadd MC.SingleTop_UL17.root *MC.ST*.root
hadd MC.TTW_UL17.root *MC.TTW*.root
hadd MC.TTZ_UL17.root *MC.TTW*.root
hadd MC.TT_UL17.root uhh2.AnalysisModuleRunner.MC.TTToHadronic_JESnominal_JERnominal_UL17.root \
    uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_JESnominal_JERnominal_UL17.root \
    uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_JESnominal_JERnominal_UL17.root
hadd MC.VV_UL17.root *MC.WW* *MC.WZ* *MC.ZZ*
hadd MC.WJets_UL17.root *MC.WJets*.root
hadd MC.ZJets_UL17.root *MC.ZJets*.root


cd $CMSSW_BASE/src/UHH2/AZH
