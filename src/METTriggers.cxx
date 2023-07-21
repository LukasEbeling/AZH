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

  HLT_PFMET120_PFMHT120_IDTight.reset(new TriggerSelection("HLT_PFMET120_PFMHT120_IDTight_v*"));
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


bool METTriggers::passes_UL16(const Event &event) {
  return HLT_PFMET120_PFMHT120_IDTight->passes(event);
}

bool METTriggers::passes_UL17(const Event &event) {
  //if (RunPeriod != "RunB") return true;
  return HLT_PFMET120_PFMHT120_IDTight->passes(event);
}

bool METTriggers::passes_UL18(const Event &event) {
  return HLT_PFMET120_PFMHT120_IDTight->passes(event);
}

bool METTriggers::passes(const Event &event) {
  bool passMET = false;

  if ((year == Year::isUL16preVFP) || (year == Year::isUL16postVFP)) {
    passMET = passes_UL16(event);
  }

  else if (year == Year::isUL17) {
    passMET = passes_UL17(event);
  }

  else if (year == Year::isUL18) {
    passMET = passes_UL18(event);
  }

  return passMET;
}