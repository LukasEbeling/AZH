#include <boost/algorithm/string.hpp>
#include <iostream>
#include <memory>
#include <string>

#include "UHH2/AZH/include/METTriggers.h"

using namespace std;
using namespace uhh2;
using namespace boost;

METTriggers::METTriggers(Context &ctx) {
  year = extract_year(ctx);

  first_event = true;

  dataset_version = ctx.get("dataset_version");
  dataset_type = ctx.get("dataset_type");

  getPeriod();

  // 2016
  HLT_PFMET300.reset(new TriggerSelection("HLT_PFMET300_v*"));
  HLT_MET200.reset(new TriggerSelection("HLT_MET200_v*"));
  HLT_PFHT300_PFMET110.reset(new TriggerSelection("HLT_PFHT300_PFMET110_v*"));
  HLT_PFMET170_HBHECleaned.reset(new TriggerSelection("HLT_PFMET170_HBHECleaned_v*"));
  HLT_PFMET120_PFMHT120_IDTight.reset(new TriggerSelection("HLT_PFMET120_PFMHT120_IDTight_v*"));
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight.reset(new TriggerSelection("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*"));
  // 2017, 2018
  HLT_PFMET200_HBHECleaned.reset(new TriggerSelection("HLT_PFMET200_HBHECleaned_v*"));
  HLT_PFMET200_HBHE_BeamHaloCleaned.reset(new TriggerSelection("HLT_PFMET200_HBHE_BeamHaloCleaned_v*"));
  HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned.reset(new TriggerSelection("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v*"));
  HLT_PFMET120_PFMHT120_IDTight_PFHT60.reset(new TriggerSelection("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*"));
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60.reset(
      new TriggerSelection("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*"));
  HLT_PFHT500_PFMET100_PFMHT100_IDTight.reset(new TriggerSelection("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*"));
  HLT_PFHT700_PFMET85_PFMHT85_IDTight.reset(new TriggerSelection("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v*"));
  HLT_PFHT800_PFMET75_PFMHT75_IDTight.reset(new TriggerSelection("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v*"));

}

void METTriggers::getPeriod() {
  RunPeriod = "";

  bool is_mc = dataset_type == "MC";

  if (!is_mc) {
    std::vector<string> fields;
    // DATA have names: PD_PERIOD_YEAR
    // Split dataset_version at "_" and return PERIOD
    split(fields, dataset_version, is_any_of("_"));
    RunPeriod = fields[1];
  }
}


//UL16
void METTriggers::passes_UL16(const Event &event, bool &passMET) {
  passMET = HLT_PFMET300->passes(event) || HLT_MET200->passes(event) || HLT_PFHT300_PFMET110->passes(event) ||
            HLT_PFMET170_HBHECleaned->passes(event) || HLT_PFMET120_PFMHT120_IDTight->passes(event) ||
            HLT_PFMETNoMu120_PFMHTNoMu120_IDTight->passes(event);
}


//UL17
void METTriggers::passes_UL17(const Event &event, bool &passMET) {
  if (RunPeriod != "RunB") {
    passMET = passMET || HLT_PFMET200_HBHECleaned->passes(event) || HLT_PFMET200_HBHE_BeamHaloCleaned->passes(event) ||
              HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned->passes(event) ||
              HLT_PFMET120_PFMHT120_IDTight_PFHT60->passes(event) ||
              HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60->passes(event);
  }
  passMET = passMET || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight->passes(event) ||
            HLT_PFHT500_PFMET100_PFMHT100_IDTight->passes(event) ||
            HLT_PFHT700_PFMET85_PFMHT85_IDTight->passes(event) || HLT_PFHT800_PFMET75_PFMHT75_IDTight->passes(event);
}


//UL18
void METTriggers::passes_UL18(const Event &event, bool &passMET) {
  passMET = HLT_PFMET200_HBHE_BeamHaloCleaned->passes(event) || 
      HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned->passes(event) ||
      HLT_PFMET120_PFMHT120_IDTight->passes(event) || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight->passes(event) ||
      HLT_PFMET120_PFMHT120_IDTight_PFHT60->passes(event) ||
      HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60->passes(event) ||
      HLT_PFHT500_PFMET100_PFMHT100_IDTight->passes(event) || HLT_PFHT700_PFMET85_PFMHT85_IDTight->passes(event) ||
      HLT_PFHT800_PFMET75_PFMHT75_IDTight->passes(event);
}

bool METTriggers::passes(const Event &event) {
  bool passMET = false;

  if ((year == Year::isUL16preVFP) || (year == Year::isUL16postVFP)) {
    passes_UL16(event, passMET);
  }

  else if (year == Year::isUL17) {
    passes_UL17(event, passMET);
  }

  else if (year == Year::isUL18) {
    passes_UL18(event, passMET);
  }

  return passMET;
}