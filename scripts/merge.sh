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

cp *AToZHToInv_MA-600_MH-400_JESnominal_JERnominal_UL17.root MC.AZH_600_400_UL17.root
cp *AToZHToInv_MA-700_MH-450_JESnominal_JERnominal_UL17.root MC.AZH_700_450_UL17.root
cp *AToZHToInv_MA-750_MH-400_JESnominal_JERnominal_UL17.root MC.AZH_750_400_UL17.root
cp *AToZHToInv_MA-750_MH-650_JESnominal_JERnominal_UL17.root MC.AZH_750_650_UL17.root
cp *AToZHToInv_MA-800_MH-400_JESnominal_JERnominal_UL17.root MC.AZH_800_400_UL17.root
cp *AToZHToInv_MA-1000_MH-400_JESnominal_JERnominal_UL17.root MC.AZH_1000_400_UL17.root
cp *AToZHToInv_MA-1000_MH-850_JESnominal_JERnominal_UL17.root MC.AZH_1000_850_UL17.root

hadd -f MC.QCD_UL17.root $(ls -S *MC.QCD*nominal*nominal*) 
hadd -f MC.SingleTop_UL17.root $(ls -S *MC.ST*nominal*nominal*)
hadd -f MC.TTW_UL17.root $(ls -S *MC.TTW*nominal*nominal*)
hadd -f MC.TTZ_UL17.root $(ls -S *MC.TTZ*nominal*nominal*)
hadd -f MC.VV_UL17.root $(ls -S *MC.WW*nominal*nominal*) $(ls -S *MC.WZ*nominal*nominal*) $(ls -S *MC.ZZ*nominal*nominal*) 
hadd -f MC.WJets_UL17.root $(ls -S *MC.WJets*nominal*nominal*)
hadd -f MC.DYJets_UL17.root $(ls -S *MC.ZJets*nominal*nominal*) $(ls -S *MC.DYJets*nominal*nominal*)
hadd -f MC.TT_UL17.root $(ls -S *MC.TTTo*nominal*nominal*)

for variation in "JESup" "JESdown" "JERup" "JERdown"
  do
    cp *AToZHToInv_MA-600_MH-400*${variation}*UL17.root MC.AZH_600_400_UL17_${variation}.root
    cp *AToZHToInv_MA-700_MH-450*${variation}*UL17.root MC.AZH_700_450_UL17_${variation}.root
    cp *AToZHToInv_MA-750_MH-400*${variation}*UL17.root MC.AZH_750_400_UL17_${variation}.root
    cp *AToZHToInv_MA-750_MH-650*${variation}*UL17.root MC.AZH_750_650_UL17_${variation}.root
    cp *AToZHToInv_MA-800_MH-400*${variation}*UL17.root MC.AZH_800_400_UL17_${variation}.root
    cp *AToZHToInv_MA-1000_MH-400*${variation}*UL17.root MC.AZH_1000_400_UL17_${variation}.root
    cp *AToZHToInv_MA-1000_MH-850*${variation}*UL17.root MC.AZH_1000_850_UL17_${variation}.root

    hadd -f MC.QCD_UL17_${variation}.root $(ls -S *MC.QCD*${variation}*) 
    hadd -f MC.SingleTop_UL17_${variation}.root $(ls -S *MC.ST*${variation}*)
    hadd -f MC.TTW_UL17_${variation}.root $(ls -S *MC.TTW*${variation}*)
    hadd -f MC.TTZ_UL17_${variation}.root $(ls -S *MC.TTZ*${variation}*)
    hadd -f MC.VV_UL17_${variation}.root $(ls -S *MC.WW*${variation}*) $(ls -S *MC.WZ*${variation}*) $(ls -S *MC.ZZ*${variation}*) 
    hadd -f MC.WJets_UL17_${variation}.root $(ls -S *MC.WJets*${variation}*)
    hadd -f MC.DYJets_UL17_${variation}.root $(ls -S *MC.ZJets*${variation}*) $(ls -S *MC.DYJets*${variation}*)
    hadd -f MC.TT_UL17_${variation}.root $(ls -S *MC.TTTo*${variation}*)
  done



cd $CMSSW_BASE/src/UHH2/AZH
