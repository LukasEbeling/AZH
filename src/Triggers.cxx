#include <boost/algorithm/string.hpp>
#include <iostream>
#include <memory>
#include <string>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/AZH/include/Utils.h"
#include "UHH2/AZH/include/Triggers.h"

using namespace std;
using namespace uhh2;
using namespace boost;


Triggers::Triggers(Context & ctx) {
  year = extract_year(ctx);

  first_event = true;

  dataset_version = ctx.get("dataset_version");
  dataset_type = ctx.get("dataset_type");

  getDataset();
  getPeriod();

  // DiEle HLTs
  HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf.reset(new TriggerSelection("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*"));
  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ.reset(new TriggerSelection("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"));
  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL.reset(new TriggerSelection("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*"));
  // DiMuon HLTs
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ.reset(new TriggerSelection("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"));
  HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ.reset(new TriggerSelection("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"));
  //Single Muon HLTs
  HLT_IsoMu24.reset(new TriggerSelection("HLT_IsoMu24_v*"));
  HLT_IsoMu24tk.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
  HLT_IsoMu24_eta2p1.reset(new TriggerSelection("HLT_IsoMu24_eta2p1_v*"));
  HLT_IsoMu22_eta2p1.reset(new TriggerSelection("HLT_IsoMu22_eta2p1_v*"));
  HLT_IsoTkMu22_eta2p1.reset(new TriggerSelection("HLT_IsoTkMu22_eta2p1_v*"));
  HLT_IsoMu27.reset(new TriggerSelection("HLT_IsoMu27_v*"));
  //Single Electron HLTs
  HLT_Ele32_eta2p1_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele32_eta2p1_WPTight_Gsf_v*"));
  HLT_Ele27_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
  HLT_Ele25_eta2p1_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele25_eta2p1_WPTight_Gsf_v*"));
  HLT_Ele32_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele32_WPTight_Gsf_v*"));
  HLT_Ele35_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele35_WPTight_Gsf_v*"));
  HLT_Ele38_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele38_WPTight_Gsf_v*"));
  HLT_Ele40_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele40_WPTight_Gsf_v*"));
  HLT_Ele32_WPTight_Gsf_L1DoubleEG.reset(new TriggerSelection("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*"));
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8.reset(new TriggerSelection("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*"));
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8.reset(new TriggerSelection("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*"));
  // Electron + Muon HLTs
  HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ.reset(new TriggerSelection("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"));
  HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL.reset(new TriggerSelection("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*"));
  HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ.reset(new TriggerSelection("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"));
  HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL.reset(new TriggerSelection("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*"));
  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ.reset(new TriggerSelection("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"));
  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL.reset(new TriggerSelection("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*"));

  // Efficiency Histograms
  h_denominator.reset(new TriggerHists(ctx, "Trigger_Denominator"));
  h_numerator.reset(new TriggerHists(ctx, "Trigger_Numerator"));
}


void Triggers::getDataset() {
  PrimaryDataset = "";

  bool is_mc = dataset_type == "MC";

  if (!is_mc){
    if (dataset_version.find("DoubleMuon") != std::string::npos) {
      PrimaryDataset = "DoubleMuon";
    } else if (dataset_version.find("DoubleEG") != std::string::npos) {
      PrimaryDataset = "DoubleEG";
    } else if (dataset_version.find("EGamma") != std::string::npos) {
      PrimaryDataset = "EGamma";
    } else if(dataset_version.find("MuonEG") != std::string::npos){
      PrimaryDataset = "MuonEG";
    } else if(dataset_version.find("SingleMuon") != std::string::npos){
      PrimaryDataset = "SingleMuon";
    } else if(dataset_version.find("SingleElectron") != std::string::npos){
      PrimaryDataset = "SingleElectron";
    } else {
      throw std::invalid_argument("Invalid PrimaryDataset. Interrupting trigger module.");
    }
  }
}


void Triggers::getPeriod() {
  RunPeriod = "";

  bool is_mc = dataset_type == "MC";

  if(!is_mc){
    std::vector<string> fields;
    // DATA have names: PD_PERIOD_YEAR
    // Split dataset_version at "_" and return PERIOD
    split(fields, dataset_version, is_any_of("_"));
    RunPeriod = fields[1];
  }
}


void Triggers::passesUL16(const Event & event, bool & passDiEle, bool & passDiMu, bool & passEleMu, bool & passSingleMu, bool & passSingleEle) {
  // Electrons
  passSingleEle = HLT_Ele32_eta2p1_WPTight_Gsf->passes(event)
              || HLT_Ele27_WPTight_Gsf->passes(event)
              || HLT_Ele25_eta2p1_WPTight_Gsf->passes(event);  // check -- TOP recommendation says last two paths might not be efficient in turn-on
  passDiEle = HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf->passes(event)
           || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ->passes(event);  // check -- TOP recommends these two

  // Muons
  passSingleMu = HLT_IsoMu24->passes(event)
              || HLT_IsoMu24tk->passes(event);  // check
  if (RunPeriod != "RunA" || RunPeriod != "RunB") {
    passSingleMu = passSingleMu
                || HLT_IsoMu22_eta2p1->passes(event)
                || HLT_IsoTkMu22_eta2p1->passes(event);  // check -- TOP lists these in two rows in a table. I assume that this signifies OR.
  }
  passDiMu = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ->passes(event)  // check
          || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ->passes(event);  // check

  // Electron / Muon
  passEleMu = HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ->passes(event);  // check -- DZ version unprescaled throughout
  bool is_mc = dataset_type == "MC";
  if (!is_mc) passEleMu = passEleMu || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ->passes(event);  // check -- DZ version unprescaled throughout; data only
  if (RunPeriod != "RunH") {
    // These two HLT paths are available only from RunH
    // cf.: https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2016
    passEleMu = passEleMu
             || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL->passes(event)  // check -- non DZ version prescaled in RunH
             || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL->passes(event);  // check -- non DZ version prescaled in RunH
  }
  if (RunPeriod == "RunF" || RunPeriod == "RunG") {
    passEleMu = passEleMu
             || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ->passes(event)
             || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL->passes(event);  // check
  }
}


void Triggers::passesUL17(const Event & event, bool & passDiEle, bool & passDiMu, bool & passEleMu, bool & passSingleMu, bool & passSingleEle) {
  // Electrons
  passSingleEle = HLT_Ele35_WPTight_Gsf->passes(event)
               || HLT_Ele38_WPTight_Gsf->passes(event)
               || HLT_Ele40_WPTight_Gsf->passes(event)
               || HLT_Ele32_WPTight_Gsf_L1DoubleEG->passes(event);  // check -- bit confusing Twiki (TOP), decided to take OR
  passDiEle = HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ->passes(event)
           || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL->passes(event);  // check -- TOP recommendations has (DZ); EG says not to use
  // Muons
  passSingleMu = HLT_IsoMu27->passes(event);
  if (RunPeriod == "RunB" || RunPeriod == "RunC" || RunPeriod == "RunD") {
    passSingleMu = passSingleMu || HLT_IsoMu24_eta2p1->passes(event);  // check -- TOP recommends OR for periods A-D
  }
  if (RunPeriod == "RunB") passDiMu = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ->passes(event);
  else {
    // HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 is available only from RunC onwards
    // cf.: https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2017
    // We don't check PD because this is true for all PDs
    // According to TOP Twiki, only HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ is available in RunB
    // cf.: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTriggerYear2017
    passDiMu = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8->passes(event)
            || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8->passes(event); // check -- TOP only recommends MassX triggers for > RunB
  }
  // Electron/Muon
  passEleMu = HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ->passes(event)
           || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL->passes(event)
           || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ->passes(event)
           || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ->passes(event);  // check
}


void Triggers::passesUL18(const Event & event, bool & passDiEle, bool & passDiMu, bool & passEleMu, bool & passSingleMu, bool & passSingleEle) {
  // Electrons
  passSingleEle = HLT_Ele32_WPTight_Gsf->passes(event);  // check -- only singleEle in 18
  passDiEle = HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL->passes(event)
           || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ->passes(event);  // check -- TOP recommendations has (DZ); EG says not to use
  // Muons
  passSingleMu = HLT_IsoMu24->passes(event);  // check -- only singleMu in 18
  passDiMu = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8->passes(event);  // check - only diMu in 18
  // Electron/Muon
  passEleMu = HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ->passes(event)
           || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL->passes(event)
           || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ->passes(event)
           || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ->passes(event); // check -- TOP recommends OR between these 4
}


bool Triggers::passes(const Event & event) {
  // Trigger list specified at https://codimd.web.cern.ch/s/WdSl0nxCa#
  bool passDiEle = false;
  bool passDiMu = false;
  bool passEleMu = false;
  bool passSingleMu = false;
  bool passSingleEle = false;

  h_denominator->fill(event);
  if((year == Year::isUL16preVFP) || (year == Year::isUL16postVFP)) {
    passesUL16(event, passDiEle, passDiMu, passEleMu, passSingleMu, passSingleEle);
  }

  else if(year == Year::isUL17) {
    passesUL17(event, passDiEle, passDiMu, passEleMu, passSingleMu, passSingleEle);
  }

  else if(year == Year::isUL18) {
    passesUL18(event, passDiEle, passDiMu, passEleMu, passSingleMu, passSingleEle);
  }

  if (passDiEle || passSingleEle || passDiMu || passSingleMu || passEleMu) h_numerator->fill(event);

  // For MC PD="", hence OR logic  on HLT is applied
  // For Data, OR logic combined with other HLT path check
  // to avoid duplicate selection of events (i.e. double counting)
  return ((PrimaryDataset=="" && (passDiEle || passDiMu || passEleMu || passSingleEle || passSingleMu)) ||
      (PrimaryDataset=="DoubleMuon" && (passDiMu)) ||
      ((PrimaryDataset=="DoubleEG" || PrimaryDataset=="EGamma") && passDiEle && !passDiMu) ||
      ((PrimaryDataset =="MuonEG") && passEleMu && !passDiEle && !passDiMu) ||
      ((PrimaryDataset =="SingleMuon") && passSingleMu && !passEleMu && !passDiEle && !passDiMu) ||
      ((PrimaryDataset =="SingleElectron") && passSingleEle && !passSingleMu && !passEleMu && !passDiEle && !passDiMu));
}

