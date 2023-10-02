#!/bin/bash

signals="/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/config/signals.txt"

cd $CMSSW_BASE/src/UHH2/AZH/data/output_02_reconstruction/MC/UL17/


for variation in "_JESup" "_JESdown" "_JERup" "_JERdown" "_JESnominal_JERnominal"; do 
    
    var_short="${variation/_JESnominal_JERnominal/}"
    
    hadd -f MC.QCD_UL17${var_short}.root $(ls -S *MC.QCD*${variation}*) 
    hadd -f MC.SingleTop_UL17${var_short}.root $(ls -S *MC.ST*${variation}*)
    hadd -f MC.TTW_UL17${var_short}.root $(ls -S *MC.TTW*${variation}*)
    hadd -f MC.TTZ_UL17${var_short}.root $(ls -S *MC.TTZ*${variation}*)
    hadd -f MC.VV_UL17${var_short}.root $(ls -S *MC.WW*${variation}*) $(ls -S *MC.WZ*${variation}*) $(ls -S *MC.ZZ*${variation}*) 
    hadd -f MC.WJets_UL17${var_short}.root $(ls -S *MC.WJets*${variation}*)
    hadd -f MC.DYJets_UL17${var_short}.root $(ls -S *MC.ZJets*${variation}*) $(ls -S *MC.DYJets*${variation}*)
    hadd -f MC.TT_UL17${var_short}.root $(ls -S *MC.TTTo*${variation}*)


    for mass in $(cat "$signals"); do
        mass_short="${mass/MA-/}"
        mass_short="${mass_short/MH-/}"
        cp *AToZHToInv_${mass}*${variation}*_UL17.root MC.AZH_${mass_short}_UL17${var_short}.root
        done

    done
