#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/AdditionalSelections.h"

#include "UHH2/AZH/include/MetFilters.h"

using namespace std;
using namespace uhh2;

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2

METFilters::METFilters(Context & ctx) {

  year = extract_year(ctx);

  // The standard primary vertex id, i.e. ndof >= 4, |z| < 24cm, rho < 2cm 
  PrimaryVertexId pvid=StandardPrimaryVertexId();
  // Apply the above primary vertex id
  modules.emplace_back(new PrimaryVertexCleaner(pvid));

  METSel.reset(new AndSelection(ctx, "metfilters"));

  if(is_UL(year)) {
    METSel->add<TriggerSelection>("goodVertices", "Flag_goodVertices");
    METSel->add<TriggerSelection>("globalSuperTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter");
    METSel->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
    METSel->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
    METSel->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
    METSel->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter");
    METSel->add<TriggerSelection>("BadPFMuonDzFilter", "Flag_BadPFMuonDzFilter");
    // Not recommended in TWiki
    // METSel->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter");
    METSel->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
    if(year == Year::isUL17 || year == Year::isUL18) {
      METSel->add<TriggerSelection>("ecalBadCalibFilter", "Flag_ecalBadCalibFilter");
      // Not recommended in TWiki
      // METSel->add<TriggerSelection>("hfNoisyHitsFilter", "Flag_hfNoisyHitsFilter");
    }
  }
  else {
    throw runtime_error("METFilters: Non-UL years not implemented");
  }

  METSel->add<NPVSelection>("1 good PV",1,-1,pvid);
}

bool METFilters::process(Event & event) {

  // Check if event passes MET filters
  if(!METSel->passes(event)) return false;
  
  // Apply cleaning of primary vertex
  // Do so iff event passes MET filters
  // Implementation inspired from CommonModules
  for(auto & m : modules){
    m->process(event);
  }

  return true;
}