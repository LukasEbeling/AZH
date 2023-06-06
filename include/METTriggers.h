#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/FlavorParticle.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/fwd.h"

#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/YearRunSwitchers.h"

#include "UHH2/AZH/include/Utils.h"

namespace uhh2 {
class METTriggers : public uhh2::Selection {
public:
  METTriggers(uhh2::Context &ctx);
  virtual bool passes(const uhh2::Event &event) override;

private:
  Year year;
  // To check if trigger exists
  bool first_event;
  std::string RunPeriod;
  std::string dataset_version;
  std::string dataset_type;

  std::unique_ptr<uhh2::Selection> HLT_PFMET120_PFMHT120_IDTight;

  // Methods
  void getDataset();
  void getPeriod();
  bool passes_UL16(const Event &event);
  bool passes_UL17(const Event &event);
  bool passes_UL18(const Event &event);
};
} // namespace uhh2