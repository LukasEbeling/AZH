#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
cd /nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src

cmsenv



cd "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/data/output_01_preselection/MC/UL17/workdir_Preselection_MC_JESnominal_JERnominal_UL17/"

# Name pattern to match
name_pattern="TTToSemi"

# Search for files containing the name pattern
found_files=$(find -type f -name "*$name_pattern*")
LD_PRELOAD=/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/startup_C.so hadd -f merged_${name_pattern}.root ${found_files[@]} 

cd "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/"
