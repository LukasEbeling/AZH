#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/FlavorParticle.h"
#include "UHH2/core/include/AnalysisModule.h"

#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/AZH/include/AtoZHHists.h"

namespace uhh2 {
  class Triggers: public uhh2::Selection {
    public:
      Triggers(uhh2::Context & ctx);
      virtual bool passes(const uhh2::Event & event) override;

    private:
      Year year;
      // To check if trigger exists
      bool first_event;
      Event::TriggerIndex trigindex_ele23, trigindex_ele27;
      Event::TriggerIndex trigindex_mudz, trigindex_mu8, trigindex_mu3p8;
      std::string PrimaryDataset;
      std::string RunPeriod;
      std::string dataset_version;
      std::string dataset_type;
      // Electrons
      std::unique_ptr<uhh2::Selection> HLT_Ele32_WPTight_Gsf;
      std::unique_ptr<uhh2::Selection> HLT_Ele27_WPTight_Gsf;
      std::unique_ptr<uhh2::Selection> HLT_Ele25_eta2p1_WPTight_Gsf;
      std::unique_ptr<uhh2::Selection> HLT_Ele32_eta2p1_WPTight_Gsf;
      std::unique_ptr<uhh2::Selection> HLT_Ele35_WPTight_Gsf;
      std::unique_ptr<uhh2::Selection> HLT_Ele38_WPTight_Gsf;
      std::unique_ptr<uhh2::Selection> HLT_Ele40_WPTight_Gsf;
      std::unique_ptr<uhh2::Selection> HLT_Ele32_WPTight_Gsf_L1DoubleEG;
      std::unique_ptr<uhh2::Selection> HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf;
      std::unique_ptr<uhh2::Selection> HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
      std::unique_ptr<uhh2::Selection> HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
      // Muons
      std::unique_ptr<uhh2::Selection> HLT_IsoMu24tk;
      std::unique_ptr<uhh2::Selection> HLT_IsoMu27;
      std::unique_ptr<uhh2::Selection> HLT_IsoMu22_eta2p1;
      std::unique_ptr<uhh2::Selection> HLT_IsoTkMu22_eta2p1;
      std::unique_ptr<uhh2::Selection> HLT_IsoMu24_eta2p1;
      std::unique_ptr<uhh2::Selection> HLT_IsoMu24;
      std::unique_ptr<uhh2::Selection> HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
      std::unique_ptr<uhh2::Selection> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, hlt_mu17tkmu8_dz_m8;
      std::unique_ptr<uhh2::Selection> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
      std::unique_ptr<uhh2::Selection> HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
      std::unique_ptr<uhh2::Selection> HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
      std::unique_ptr<uhh2::Selection> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
      // Electron / Muon
      std::unique_ptr<uhh2::Selection> HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
      std::unique_ptr<uhh2::Selection> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
      std::unique_ptr<uhh2::Selection> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
      std::unique_ptr<uhh2::Selection> HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
      void getDataset();
      void getPeriod();
      void passesUL16(const Event & event, bool & passDiEle,  bool & passDiMu, bool & passEleMu, bool & passSingleMu, bool & passSingleEle);
      void passesUL17(const Event & event, bool & passDiEle,  bool & passDiMu, bool & passEleMu, bool & passSingleMu, bool & passSingleEle);
      void passesUL18(const Event & event, bool & passDiEle,  bool & passDiMu, bool & passEleMu, bool & passSingleMu, bool & passSingleEle);

      // Trigger Histograms
      std::unique_ptr<Hists> h_denominator;
      std::unique_ptr<Hists> h_numerator;
    };
}
