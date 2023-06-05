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

  // Triggers
  // 2016
  std::unique_ptr<uhh2::Selection> HLT_PFMET300;
  std::unique_ptr<uhh2::Selection> HLT_MET200;
  std::unique_ptr<uhh2::Selection> HLT_PFHT300_PFMET110;
  std::unique_ptr<uhh2::Selection> HLT_PFMET170_HBHECleaned;
  std::unique_ptr<uhh2::Selection> HLT_PFMET120_PFMHT120_IDTight; // also 17, 18
  std::unique_ptr<uhh2::Selection> HLT_PFMET120NoMu_PFMHT120NoMu_IDTight;
  // 2017, 2018
  std::unique_ptr<uhh2::Selection> HLT_PFMET200_HBHECleaned;                 // MET POG
  std::unique_ptr<uhh2::Selection> HLT_PFMET200_HBHE_BeamHaloCleaned;        // MET POG
  std::unique_ptr<uhh2::Selection> HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned; // MET POG
  std::unique_ptr<uhh2::Selection> HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
  std::unique_ptr<uhh2::Selection> HLT_PFMET120_PFMHT120_IDTight_PFHT60;
  std::unique_ptr<uhh2::Selection> HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  std::unique_ptr<uhh2::Selection> HLT_PFHT500_PFMET100_PFMHT100_IDTight;
  std::unique_ptr<uhh2::Selection> HLT_PFHT700_PFMET85_PFMHT85_IDTight;
  std::unique_ptr<uhh2::Selection> HLT_PFHT800_PFMET75_PFMHT75_IDTight;

  // Methods
  void getDataset();
  void getPeriod();
  void passes_UL16(const Event &event, bool &passMET);
  void passes_UL17(const Event &event, bool &passMET);
  void passes_UL18(const Event &event, bool &passMET);
};
} // namespace uhh2