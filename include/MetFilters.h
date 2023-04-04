#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
class METFilters: public uhh2::AnalysisModule {
public:
  explicit METFilters(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event);

private:
  Year year;
  std::vector<std::unique_ptr<AnalysisModule>> modules;
  std::unique_ptr<uhh2::AndSelection> METSel;
};