#include <iostream>
#include <memory>

#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/AZH/include/ScaleFactors.h"
#include "UHH2/AZH/include/Utils.h"

#include <TFile.h>
#include <TH2D.h>
#include <TH2F.h>

using namespace std;
using namespace uhh2;

//-----------------------------------------------------------------------
// Leptons Scale Factors
//-----------------------------------------------------------------------
// SFs list and files taken from:
// https://github.com/UHH2/UHH2/blob/RunII_106X_v2/common/src/LeptonScaleFactors.cxx

const string uhh_data_path = string(std::getenv("CMSSW_BASE")) + "/src/UHH2/common/UHH2-data/";

ScaleFactors2016preVFP::ScaleFactors2016preVFP(Context &ctx) {

  string filepath_sf_mu_reco = (string)getenv("CMSSW_BASE") +
                                "/src/UHH2/AZH/data/ScaleFactors/reco/Muon_Reconstruction_SF_UL16.root";

  mu_sf_reco.reset(new MCMuonScaleFactor(ctx, filepath_sf_mu_reco, "Reconstruction_SF", 0.0,
                                         "reco", false, "nominal", "muons", true));

  string filepath_sf_mu_id =
      uhh_data_path + "/muon_SFs_UL/UL16preVFP/" + "Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root";

  mu_sf_id.reset(new MCMuonScaleFactor(ctx, filepath_sf_mu_id, "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0.0,
                                       "tight_id", false, "nominal", "muons", true));

  string filepath_sf_iso =
      uhh_data_path + "/muon_SFs_UL/UL16preVFP/" + "Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root";

  mu_sf_iso.reset(new MCMuonScaleFactor(ctx, filepath_sf_iso, "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", 0.0,
                                        "isolation", false, "nominal", "muons", true));

  string filepath_sf_ele_id =
      uhh_data_path + "/egamma_SFs_UL/UL16preVFP/" + "egammaEffi.txt_Ele_wp80iso_preVFP_EGM2D.root";

  ele_sf_id.reset(
      new MCElecScaleFactor(ctx, filepath_sf_ele_id, 0.0, "tight_id", "nominal", "electrons", "EGamma_SF2D", false));

  string filepath_sf_reco =
      uhh_data_path + "/egamma_SFs_UL/UL16preVFP/" + "egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root";

  ele_sf_reco.reset(
      new MCElecScaleFactor(ctx, filepath_sf_reco, 0.0, "reco", "nominal", "electrons", "EGamma_SF2D", false));
}

bool ScaleFactors2016preVFP::process(Event &event) {

  mu_sf_reco->process(event);
  mu_sf_id->process(event);
  mu_sf_iso->process(event);

  ele_sf_id->process(event);
  ele_sf_reco->process(event);

  return true;
}

ScaleFactors2016postVFP::ScaleFactors2016postVFP(Context &ctx) {

  string filepath_sf_mu_reco = (string)getenv("CMSSW_BASE") +
                                "/src/UHH2/AZH/data/ScaleFactors/reco/Muon_Reconstruction_SF_UL16.root";

  mu_sf_reco.reset(new MCMuonScaleFactor(ctx, filepath_sf_mu_reco, "Reconstruction_SF", 0.0,
                                         "reco", false, "nominal", "muons", true));

  string filepath_sf_mu_id =
      uhh_data_path + "/muon_SFs_UL/UL16postVFP/" + "Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root";

  mu_sf_id.reset(new MCMuonScaleFactor(ctx, filepath_sf_mu_id, "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0.0,
                                       "tight_id", false, "nominal", "muons", true));

  string filepath_sf_iso =
      uhh_data_path + "/muon_SFs_UL/UL16postVFP/" + "Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root";

  mu_sf_iso.reset(new MCMuonScaleFactor(ctx, filepath_sf_iso, "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", 0.0,
                                        "isolation", false, "nominal", "muons", true));

  string filepath_sf_ele_id =
      uhh_data_path + "/egamma_SFs_UL/UL16postVFP/" + "egammaEffi.txt_Ele_wp80iso_postVFP_EGM2D.root";

  ele_sf_id.reset(
      new MCElecScaleFactor(ctx, filepath_sf_ele_id, 0.0, "tight_id", "nominal", "electrons", "EGamma_SF2D", false));

  string filepath_sf_reco =
      uhh_data_path + "/egamma_SFs_UL/UL16postVFP/" + "egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root";

  ele_sf_reco.reset(
      new MCElecScaleFactor(ctx, filepath_sf_reco, 0.0, "reco", "nominal", "electrons", "EGamma_SF2D", false));
}

bool ScaleFactors2016postVFP::process(Event &event) {

  mu_sf_reco->process(event);
  mu_sf_id->process(event);
  mu_sf_iso->process(event);

  ele_sf_id->process(event);
  ele_sf_reco->process(event);

  return true;
}

ScaleFactors2017::ScaleFactors2017(Context &ctx) {

  string filepath_sf_mu_reco = (string)getenv("CMSSW_BASE") +
                                "/src/UHH2/AZH/data/ScaleFactors/reco/Muon_Reconstruction_SF_UL17.root";

  mu_sf_reco.reset(new MCMuonScaleFactor(ctx, filepath_sf_mu_reco, "Reconstruction_SF", 0.0,
                                         "reco", false, "nominal", "muons", true));

  string filepath_sf_mu_id =
      uhh_data_path + "/muon_SFs_UL/UL17/" + "Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root";

  mu_sf_id.reset(new MCMuonScaleFactor(ctx, filepath_sf_mu_id, "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0.0,
                                       "tight_id", false, "nominal", "muons", true));

  string filepath_sf_iso =
      uhh_data_path + "/muon_SFs_UL/UL17/" + "Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root";

  mu_sf_iso.reset(new MCMuonScaleFactor(ctx, filepath_sf_iso, "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", 0.0,
                                        "isolation", false, "nominal", "muons", true));

  string filepath_sf_ele_id = uhh_data_path + "/egamma_SFs_UL/UL17/" + "egammaEffi.txt_EGM2D_MVA80iso_UL17.root";

  ele_sf_id.reset(
      new MCElecScaleFactor(ctx, filepath_sf_ele_id, 0.0, "tight_id", "nominal", "electrons", "EGamma_SF2D", false));

  string filepath_sf_reco = uhh_data_path + "/egamma_SFs_UL/UL17/" + "egammaEffi_ptAbove20.txt_EGM2D_UL2017.root";

  ele_sf_reco.reset(
      new MCElecScaleFactor(ctx, filepath_sf_reco, 0.0, "reco", "nominal", "electrons", "EGamma_SF2D", false));
}

bool ScaleFactors2017::process(Event &event) {

  mu_sf_reco->process(event);
  mu_sf_id->process(event);
  mu_sf_iso->process(event);

  ele_sf_id->process(event);
  ele_sf_reco->process(event);

  return true;
}

ScaleFactors2018::ScaleFactors2018(Context &ctx) {

  string filepath_sf_mu_reco = (string)getenv("CMSSW_BASE") +
                                "/src/UHH2/AZH/data/ScaleFactors/reco/Muon_Reconstruction_SF_UL18.root";

  mu_sf_reco.reset(new MCMuonScaleFactor(ctx, filepath_sf_mu_reco, "Reconstruction_SF", 0.0,
                                         "reco", false, "nominal", "muons", true));

  string filepath_sf_mu_id =
      uhh_data_path + "/muon_SFs_UL/UL18/" + "Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root";

  mu_sf_id.reset(new MCMuonScaleFactor(ctx, filepath_sf_mu_id, "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0.0,
                                       "tight_id", false, "nominal", "muons", true));

  string filepath_sf_iso =
      uhh_data_path + "/muon_SFs_UL/UL18/" + "Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root";

  mu_sf_iso.reset(new MCMuonScaleFactor(ctx, filepath_sf_iso, "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", 0.0,
                                        "isolation", false, "nominal", "muons", true));

  string filepath_sf_ele_id = uhh_data_path + "/egamma_SFs_UL/UL18/" + "egammaEffi.txt_Ele_wp80iso_EGM2D.root";

  ele_sf_id.reset(
      new MCElecScaleFactor(ctx, filepath_sf_ele_id, 0.0, "tight_id", "nominal", "electrons", "EGamma_SF2D", false));

  string filepath_sf_reco = uhh_data_path + "/egamma_SFs_UL/UL18/" + "egammaEffi_ptAbove20.txt_EGM2D_UL2018.root";

  ele_sf_reco.reset(
      new MCElecScaleFactor(ctx, filepath_sf_reco, 0.0, "reco", "nominal", "electrons", "EGamma_SF2D", false));
}

bool ScaleFactors2018::process(Event &event) {

  mu_sf_reco->process(event);
  mu_sf_id->process(event);
  mu_sf_iso->process(event);

  ele_sf_id->process(event);
  ele_sf_reco->process(event);

  return true;
}

LeptonScaleFactors::LeptonScaleFactors(Context &ctx) {

  m_sf_lepton.reset(new YearSwitcher(ctx));

  m_sf_lepton->setupUL16preVFP(std::make_shared<ScaleFactors2016preVFP>(ctx));
  m_sf_lepton->setupUL16postVFP(std::make_shared<ScaleFactors2016postVFP>(ctx));
  m_sf_lepton->setupUL17(std::make_shared<ScaleFactors2017>(ctx));
  m_sf_lepton->setupUL18(std::make_shared<ScaleFactors2018>(ctx));
}

bool LeptonScaleFactors::process(Event &event) {

  m_sf_lepton->process(event);

  return true;
}

//-----------------------------------------------------------------------
// L1 prefiring weight correction
//-----------------------------------------------------------------------

L1PrefiringWeight::L1PrefiringWeight(Context &ctx) {

  year = extract_year(ctx);
  string syst_direction_ = ctx.get("SystDirection_Prefiring", "nominal");

  if (syst_direction_ == "up") {
    syst_direction = 1;
  } else if (syst_direction_ == "down") {
    syst_direction = -1;
  } else {
    syst_direction = 0;
  }
}

bool L1PrefiringWeight::process(Event &event) {

  // In UL we have L1 prefiring SF for both electrons and muons
  // In 2018 there is only muons weight
  // See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe
  // TODO: Maybe worth saving handles for all these
  if (!event.isRealData) {
    if (syst_direction == 1) {
      event.weight *= event.prefiringWeightUp;
    } else if (syst_direction == -1) {
      event.weight *= event.prefiringWeightDown;
    } else {
      event.weight *= event.prefiringWeight;
    }
  }

  return true;
}

//-------
// PUID SFs
// -------
PUIDScaleFactors::PUIDScaleFactors(Context &ctx) {
  year = extract_year(ctx);
  is_mc = ctx.get("dataset_type") == "MC";
  string filename = (string)getenv("CMSSW_BASE") +
                    "/src/UHH2/AZH/data/ScaleFactors/PUID_106XTraining_ULRun2_EffSFandUncties_v1.root";
  file_sf.reset(new TFile(filename.c_str(), "READ"));
  handle_sf_puid = ctx.declare_event_output<float>("weight_puid_sf");
  handle_sf_puid_up = ctx.declare_event_output<float>("weight_puid_sf_up");
  handle_sf_puid_down = ctx.declare_event_output<float>("weight_puid_sf_down");
}

double PUIDScaleFactors::load_sf(const Jet &jet, string variation) {
  // open root file and return value from correct histogram and bin
  std::string th2f_name = "h2_eff_sfUL";
  switch (year) {
  default:
    break;
  case Year::isUL16preVFP:
    th2f_name += "2016_T";
    break;
  case Year::isUL16postVFP:
    th2f_name += "2016APV_T";
    break;
  case Year::isUL17:
    th2f_name += "2017_T";
    break;
  case Year::isUL18:
    th2f_name += "2018_T";
    break;
  }

  TH2D *histo_sf = (TH2D *)file_sf->Get(th2f_name.c_str());
  TH2D *histo_sf_sys = (TH2D *)file_sf->Get((th2f_name + "_Systuncty").c_str());

  double pt = jet.pt();
  double eta = jet.eta();
  int bin_lead = histo_sf->GetXaxis()->FindBin(pt);
  int bin_sublead = histo_sf->GetYaxis()->FindBin(eta);

  // Define scale factors and corresponding errors
  if (variation == "nominal")
    return (double)histo_sf->GetBinContent(bin_lead, bin_sublead);
  if (variation == "up") {
    double nominal = (double)histo_sf->GetBinContent(bin_lead, bin_sublead);
    double variation = (double)histo_sf_sys->GetBinContent(bin_lead, bin_sublead);
    return nominal + variation;
  }
  if (variation == "down") {
    double nominal = (double)histo_sf->GetBinContent(bin_lead, bin_sublead);
    double variation = (double)histo_sf_sys->GetBinContent(bin_lead, bin_sublead);
    return nominal - variation;
  }
  return 1;
}

bool PUIDScaleFactors::process(Event &event) {
  double sf_nominal = 1;
  double sf_up = 1;
  double sf_down = 1;

  if (is_mc) {
    // iterate over jet collection
    for (unsigned i = 0; i < event.jets->size(); i++) {
      // skip jets with pt>50GeV
      if (event.jets->at(i).pt() >= 50)
        continue;
      // skip jets without a match to a genjet
      for (unsigned j = 0; j < event.genjets->size(); j++) {
        double dR = deltaR(event.jets->at(i).v4(), event.genjets->at(j).v4());
        if (dR < 0.4) {
          sf_nominal *= load_sf(event.jets->at(i), "nominal");
          sf_up *= load_sf(event.jets->at(i), "up");
          sf_down *= load_sf(event.jets->at(i), "down");
          break;
        }
      }
    }
    // Apply nominal SF to event.weight
    event.weight *= sf_nominal;
  }

  // Set event handles
  event.set(handle_sf_puid, sf_nominal);
  event.set(handle_sf_puid_up, sf_up);
  event.set(handle_sf_puid_down, sf_down);

  return true;
}

//-----------------------------------------------------------------------
// Trigger Scale Factors
//-----------------------------------------------------------------------
// Taken for dilepton triggers from:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Dilepton_triggers

DiEleTriggerScaleFactors::DiEleTriggerScaleFactors(Context &ctx) {

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  h_ee_trigger_sf = ctx.get_handle<float>("weight_trigger_sf_ee");
  h_ee_trigger_sf_up = ctx.get_handle<float>("weight_trigger_sf_ee_up");
  h_ee_trigger_sf_dn = ctx.get_handle<float>("weight_trigger_sf_ee_dn");

  if (!is_mc) {
    cout << "Warning: DiEleTriggerScaleFactors will not have an effect on this "
            "non-MC sample (dataset_type = '" +
                dataset_type + "')"
         << endl;
    return;
  }

  year = extract_year(ctx);

  string syst_direction_ = ctx.get("SystDirection_TriggerSF", "nominal");

  if (syst_direction_ == "up") {
    syst_direction = 1;
  } else if (syst_direction_ == "down") {
    syst_direction = -1;
  } else {
    syst_direction = 0;
  }
}

bool DiEleTriggerScaleFactors::process(Event &event) {

  if (event.isRealData) {
    event.set(h_ee_trigger_sf, 1.);
    event.set(h_ee_trigger_sf_dn, 1.);
    event.set(h_ee_trigger_sf_up, 1.);
    return true;
  }

  string filename = (string)getenv("CMSSW_BASE") + "/src/UHH2/AZH/data/ScaleFactors/triggers/";
  switch (year) {
  default:
    break;
  case Year::isUL16preVFP:
    filename += "TriggerSF_Ele_UL16.root";
    break;
  case Year::isUL16postVFP:
    filename += "TriggerSF_Ele_UL16.root";
    break;
  case Year::isUL17:
    filename += "TriggerSF_Ele_UL17.root";
    break;
  case Year::isUL18:
    filename += "TriggerSF_Ele_UL18.root";
    break;
  }

  std::unique_ptr<TFile> file_sf;
  file_sf.reset(new TFile(filename.c_str(), "READ"));

  TH2D *histo_ee_sf = (TH2D *)file_sf->Get("h2D_SF_central");
  TH2D *histo_ee_sf_up = (TH2D *)file_sf->Get("h2D_SF_up");
  TH2D *histo_ee_sf_down = (TH2D *)file_sf->Get("h2D_SF_down");

  // Assuming particle collections are pT ordered
  const unsigned leading_idx = (event.electrons->at(0).pt() > event.electrons->at(1).pt()) ? 0 : 1;
  const unsigned subleading_idx = (leading_idx == 0) ? 1 : 0;
  const Electron ele_lead = event.electrons->at(leading_idx);
  const Electron ele_sublead = event.electrons->at(subleading_idx);

  double pT_lead = ele_lead.pt();
  double pT_sublead = ele_sublead.pt();

  // Define scale factors and corresponding errors
  double sf, sf_up, sf_down;

  // Get bins corresponding to pT_ele1 and pT_ele2
  int bin_lead = histo_ee_sf->GetXaxis()->FindBin(pT_lead);
  int bin_sublead = histo_ee_sf->GetYaxis()->FindBin(pT_sublead);

  // Get corresponding SF value
  sf = histo_ee_sf->GetBinContent(bin_lead, bin_sublead);
  // Get up/dn stat errors
  sf_up = histo_ee_sf_up->GetBinContent(bin_lead, bin_sublead);
  sf_down = histo_ee_sf_down->GetBinContent(bin_lead, bin_sublead);

  event.set(h_ee_trigger_sf, sf);
  event.set(h_ee_trigger_sf_up, sf_up);
  event.set(h_ee_trigger_sf_dn, sf_down);

  if (syst_direction == 1) {
    event.weight *= sf_up;
  } else if (syst_direction == -1) {
    event.weight *= sf_down;
  } else {
    event.weight *= sf;
  }

  return true;
}

DiMuTriggerScaleFactors::DiMuTriggerScaleFactors(Context &ctx) {

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  h_mumu_trigger_sf = ctx.get_handle<float>("weight_trigger_sf_mumu");
  h_mumu_trigger_sf_up = ctx.get_handle<float>("weight_trigger_sf_mumu_up");
  h_mumu_trigger_sf_dn = ctx.get_handle<float>("weight_trigger_sf_mumu_dn");

  if (!is_mc) {
    cout << "Warning: DiMuTriggerScaleFactors will not have an effect on this "
            "non-MC sample (dataset_type = '" +
                dataset_type + "')"
         << endl;
    return;
  }

  year = extract_year(ctx);

  // Use same syst_direction as for DiEle
  // These are the same SF but for different channels
  // hence we want consistency
  string syst_direction_ = ctx.get("SystDirection_TriggerSF", "nominal");
  if (syst_direction_ == "up") {
    syst_direction = 1;
  } else if (syst_direction_ == "down") {
    syst_direction = -1;
  } else {
    syst_direction = 0;
  }
}

bool DiMuTriggerScaleFactors::process(Event &event) {

  if (event.isRealData) {
    event.set(h_mumu_trigger_sf, 1.);
    event.set(h_mumu_trigger_sf_dn, 1.);
    event.set(h_mumu_trigger_sf_up, 1.);
    return true;
  }

  string filename = (string)getenv("CMSSW_BASE") + "/src/UHH2/AZH/data/ScaleFactors/triggers/";
  switch (year) {
  default:
    break;
  case Year::isUL16preVFP:
    filename += "TriggerSF_Mu_UL16.root";
    break;
  case Year::isUL16postVFP:
    filename += "TriggerSF_Mu_UL16.root";
    break;
  case Year::isUL17:
    filename += "TriggerSF_Mu_UL17.root";
    break;
  case Year::isUL18:
    filename += "TriggerSF_Mu_UL18.root";
    break;
  }

  std::unique_ptr<TFile> file_sf;
  file_sf.reset(new TFile(filename.c_str(), "READ"));

  TH2D *histo_mumu_sf = (TH2D *)file_sf->Get("h2D_SF_central");
  TH2D *histo_mumu_sf_up = (TH2D *)file_sf->Get("h2D_SF_up");
  TH2D *histo_mumu_sf_down = (TH2D *)file_sf->Get("h2D_SF_down");

  // Assuming particle collections are pT ordered
  const unsigned leading_idx = (event.muons->at(0).pt() > event.muons->at(1).pt()) ? 0 : 1;
  const unsigned subleading_idx = (leading_idx == 0) ? 1 : 0;
  const Muon mu_lead = event.muons->at(leading_idx);
  const Muon mu_sublead = event.muons->at(subleading_idx);

  double pT_lead = mu_lead.pt();
  double pT_sublead = mu_sublead.pt();

  // Define scale factors and corresponding errors
  double sf, sf_up, sf_down;

  // Get bins corresponding to pT_ele1 and pT_ele2
  int bin_lead = histo_mumu_sf->GetXaxis()->FindBin(pT_lead);
  int bin_sublead = histo_mumu_sf->GetYaxis()->FindBin(pT_sublead);

  // Get corresponding SF value
  sf = histo_mumu_sf->GetBinContent(bin_lead, bin_sublead);
  sf_up = histo_mumu_sf_up->GetBinContent(bin_lead, bin_sublead);
  sf_down = histo_mumu_sf_down->GetBinContent(bin_lead, bin_sublead);

  event.set(h_mumu_trigger_sf, sf);
  event.set(h_mumu_trigger_sf_up, sf_up);
  event.set(h_mumu_trigger_sf_dn, sf_down);

  if (syst_direction == 1) {
    event.weight *= sf_up;
  } else if (syst_direction == -1) {
    event.weight *= sf_down;
  } else {
    event.weight *= sf;
  }

  return true;
}

EleMuTriggerScaleFactors::EleMuTriggerScaleFactors(Context &ctx) {

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  h_emu_trigger_sf = ctx.get_handle<float>("weight_trigger_sf_emu");
  h_emu_trigger_sf_up = ctx.get_handle<float>("weight_trigger_sf_emu_up");
  h_emu_trigger_sf_dn = ctx.get_handle<float>("weight_trigger_sf_emu_dn");

  if (!is_mc) {
    cout << "Warning: EleMuTriggerScaleFactors will not have an effect on this "
            "non-MC sample (dataset_type = '" +
                dataset_type + "')"
         << endl;
    return;
  }

  year = extract_year(ctx);

  // Use same syst_direction as for DiEle
  // These are the same SF but for different channels
  // hence we want consistency
  string syst_direction_ = ctx.get("SystDirection_TriggerSF", "nominal");
  if (syst_direction_ == "up") {
    syst_direction = 1;
  } else if (syst_direction_ == "down") {
    syst_direction = -1;
  } else {
    syst_direction = 0;
  }
}

bool EleMuTriggerScaleFactors::process(Event &event) {

  if (event.isRealData) {
    event.set(h_emu_trigger_sf, 1.);
    event.set(h_emu_trigger_sf_dn, 1.);
    event.set(h_emu_trigger_sf_up, 1.);
    return true;
  }

  string filename = (string)getenv("CMSSW_BASE") + "/src/UHH2/AZH/data/ScaleFactors/triggers/";
  switch (year) {
  default:
    break;
  case Year::isUL16preVFP:
    filename += "TriggerSF_2016preVFP_UL.root";
    break;
  case Year::isUL16postVFP:
    filename += "TriggerSF_2016postVFP_UL.root";
    break;
  case Year::isUL17:
    filename += "TriggerSF_2017_UL.root";
    break;
  case Year::isUL18:
    filename += "TriggerSF_2018_UL.root";
    break;
  }

  std::unique_ptr<TFile> file_sf;
  file_sf.reset(new TFile(filename.c_str(), "READ"));

  TH2F *histo_emu_sf = (TH2F *)file_sf->Get("h2D_SF_emu_lepABpt_FullError");

  // Take leading electron and leading muon
  // TODO: Cross check.
  const Muon mu_lead = event.muons->at(0);
  const Electron ele_lead = event.electrons->at(0);

  double pT_lead = -999;
  double pT_sublead = -999;

  if (mu_lead.pt() > ele_lead.pt()) {
    pT_lead = mu_lead.pt();
    pT_sublead = ele_lead.pt();
  } else {
    pT_sublead = mu_lead.pt();
    pT_lead = ele_lead.pt();
  }

  // Define scale factors and corresponding errors
  // TP is set to 0.0, but could be used for a 2% sys on T&P
  // Similarly to other methods the MCWeight of UHH2
  double sf, sf_up, sf_down;
  double stat_up = -1., stat_down = -1., tp = 0.0;
  double total_up = -1., total_down = -1.;

  // Get bins corresponding to pT_ele1 and pT_ele2
  int bin_lead = histo_emu_sf->GetXaxis()->FindBin(pT_lead);
  int bin_sublead = histo_emu_sf->GetYaxis()->FindBin(pT_sublead);

  // Get corresponding SF value
  sf = histo_emu_sf->GetBinContent(bin_lead, bin_sublead);
  // Get up/dn stat errors (here symmetric)
  stat_up = histo_emu_sf->GetBinErrorUp(bin_lead, bin_sublead);
  stat_down = histo_emu_sf->GetBinErrorLow(bin_lead, bin_sublead);
  // Total error and up/dn SFs
  total_up = sqrt(pow(stat_up, 2) + pow(tp, 2));
  total_down = sqrt(pow(stat_down, 2) + pow(tp, 2));
  sf_up = sf + total_up;
  sf_down = sf - total_down;

  event.set(h_emu_trigger_sf, sf);
  event.set(h_emu_trigger_sf_up, sf_up);
  event.set(h_emu_trigger_sf_dn, sf_down);

  if (syst_direction == 1) {
    event.weight *= sf_up;
  } else if (syst_direction == -1) {
    event.weight *= sf_down;
  } else {
    event.weight *= sf;
  }

  return true;
}

InitTriggerScaleFactors::InitTriggerScaleFactors(Context &ctx) {

  handle_channel = ctx.get_handle<int>("channel");

  m_handles.push_back(ctx.declare_event_output<float>("weight_trigger_sf_ee"));
  m_handles.push_back(ctx.declare_event_output<float>("weight_trigger_sf_ee_up"));
  m_handles.push_back(ctx.declare_event_output<float>("weight_trigger_sf_ee_dn"));
  m_handles.push_back(ctx.declare_event_output<float>("weight_trigger_sf_mumu"));
  m_handles.push_back(ctx.declare_event_output<float>("weight_trigger_sf_mumu_up"));
  m_handles.push_back(ctx.declare_event_output<float>("weight_trigger_sf_mumu_dn"));
  m_handles.push_back(ctx.declare_event_output<float>("weight_trigger_sf_emu"));
  m_handles.push_back(ctx.declare_event_output<float>("weight_trigger_sf_emu_up"));
  m_handles.push_back(ctx.declare_event_output<float>("weight_trigger_sf_emu_dn"));
}

bool InitTriggerScaleFactors::process(Event &event) {

  for (uint i = 0; i < m_handles.size(); i++) {
    event.set(m_handles.at(i), -1.);
  }

  return true;
}

TriggerScaleFactors::TriggerScaleFactors(Context &ctx) {

  handle_channel = ctx.get_handle<int>("channel");

  m_sf_dummy.reset(new InitTriggerScaleFactors(ctx));
  m_sf_ee_trigger.reset(new DiEleTriggerScaleFactors(ctx));
  m_sf_mumu_trigger.reset(new DiMuTriggerScaleFactors(ctx));
  m_sf_emu_trigger.reset(new EleMuTriggerScaleFactors(ctx));
}

bool TriggerScaleFactors::process(Event &event) {

  Channel channel = (Channel)event.get(handle_channel);

  m_sf_dummy->process(event);

  if (channel == Channel::diElectron) {
    m_sf_ee_trigger->process(event);
  }
  if (channel == Channel::diMuon) {
    m_sf_mumu_trigger->process(event);
  }
  if (channel == Channel::ElectronMuon) {
    m_sf_emu_trigger->process(event);
  }

  return true;
}

//-----------------------------------------------------------------------
// m(top) AND pT(V) CORRECTIONS
//-----------------------------------------------------------------------
//________________________________________________________________________
// https://github.com/MatthiesC/LegacyTopTagging/blob/master/src/Utils.cxx

TopPtReweighting::TopPtReweighting(Context &ctx, const bool apply) : fApply(apply) {
  const string config = ctx.get("SystDirection_TopPt", "nominal");
  h_weight_nominal = ctx.declare_event_output<float>("weight_toppt");
  h_weight_a_up = ctx.declare_event_output<float>("weight_toppt_a_up");
  h_weight_a_down = ctx.declare_event_output<float>("weight_toppt_a_down");
  h_weight_b_up = ctx.declare_event_output<float>("weight_toppt_b_up");
  h_weight_b_down = ctx.declare_event_output<float>("weight_toppt_b_down");
  h_weight_applied = ctx.declare_event_output<float>("weight_toppt_applied");
  proc_is_TT = (ctx.get("dataset_version").find("TTTo") == 0);
  if (config == "nominal") {
    applied_variation = TopPtVariation::nominal;
  } else if (config == "a_up") {
    applied_variation = TopPtVariation::a_up;
  } else if (config == "a_down") {
    applied_variation = TopPtVariation::a_down;
  } else if (config == "b_up") {
    applied_variation = TopPtVariation::b_up;
  } else if (config == "b_down") {
    applied_variation = TopPtVariation::b_down;
  } else {
    throw invalid_argument("TopPtReweighting: Invalid systematic variation given in XML config.");
  }
}

void TopPtReweighting::set_dummy_weights(Event &event) {
  event.set(h_weight_nominal, fDummyWeight);
  event.set(h_weight_a_up, fDummyWeight);
  event.set(h_weight_a_down, fDummyWeight);
  event.set(h_weight_b_up, fDummyWeight);
  event.set(h_weight_b_down, fDummyWeight);
  event.set(h_weight_applied, fDummyWeight);
}

bool TopPtReweighting::process(Event &event) {
  if (event.isRealData) {
    set_dummy_weights(event);
    return true;
  }

  unsigned int n_tops(0);
  float top_pt(-1.);
  float antitop_pt(-1.);
  for (const GenParticle &gp : *event.genparticles) {
    if (abs(gp.pdgId()) == 6)
      ++n_tops;
    if (gp.pdgId() == 6) {
      top_pt = gp.v4().pt();
    } else if (gp.pdgId() == -6) {
      antitop_pt = gp.v4().pt();
    }
  }
  if (fPtCutOff_b) {
    top_pt = min(top_pt, fPtCutOff);
    antitop_pt = min(antitop_pt, fPtCutOff);
  }

  if (!proc_is_TT) {
    set_dummy_weights(event);
    return true;
  }
  if (top_pt < 0 || antitop_pt < 0)
    throw runtime_error("TopPtReweighting::process(): top or antitop pT is negative.");

  const float sf_top_nominal = exp(fA + fB * top_pt);
  const float sf_top_a_up = exp(fA_up + fB * top_pt);
  const float sf_top_a_down = exp(fA_down + fB * top_pt);
  const float sf_top_b_up = exp(fA + fB_up * top_pt);
  const float sf_top_b_down = exp(fA + fB_down * top_pt);
  const float sf_antitop_nominal = exp(fA + fB * antitop_pt);
  const float sf_antitop_a_up = exp(fA_up + fB * antitop_pt);
  const float sf_antitop_a_down = exp(fA_down + fB * antitop_pt);
  const float sf_antitop_b_up = exp(fA + fB_up * antitop_pt);
  const float sf_antitop_b_down = exp(fA + fB_down * antitop_pt);

  const float sf_nominal = sqrt(sf_top_nominal * sf_antitop_nominal);
  const float sf_a_up = sqrt(sf_top_a_up * sf_antitop_a_up);
  const float sf_a_down = sqrt(sf_top_a_down * sf_antitop_a_down);
  const float sf_b_up = sqrt(sf_top_b_up * sf_antitop_b_up);
  const float sf_b_down = sqrt(sf_top_b_down * sf_antitop_b_down);
  float sf_applied = 1.0f;

  if (fApply) {
    if (applied_variation == TopPtVariation::nominal) {
      sf_applied *= sf_nominal;
    } else if (applied_variation == TopPtVariation::a_up) {
      sf_applied *= sf_a_up;
    } else if (applied_variation == TopPtVariation::a_down) {
      sf_applied *= sf_a_down;
    } else if (applied_variation == TopPtVariation::b_up) {
      sf_applied *= sf_b_up;
    } else if (applied_variation == TopPtVariation::b_down) {
      sf_applied *= sf_b_down;
    }
  }
  event.weight *= sf_applied;

  event.set(h_weight_nominal, sf_nominal);
  event.set(h_weight_a_up, sf_a_up);
  event.set(h_weight_a_down, sf_a_down);
  event.set(h_weight_b_up, sf_b_up);
  event.set(h_weight_b_down, sf_b_down);
  event.set(h_weight_applied, sf_applied);

  return true;
}

//___________________________________________________________________________________________
// Adapted from
// https://github.com/MatthiesC/HighPtSingleTop/blob/master/src/TheoryCorrections.cxx

VJetsReweighting::VJetsReweighting(Context &ctx, const string &weight_name)
    : is_WJets(ctx.get("dataset_version").find("WJets") == 0),
      is_ZJets(ctx.get("dataset_version").find("ZJets") == 0),
      is_DYJets(ctx.get("dataset_version").find("DYJets") == 0),
      apply_EWK(string2bool(ctx.get("VJetsReweighting_do_EWK"))),
      apply_QCD_EWK(string2bool(ctx.get("VJetsReweighting_do_QCD_EWK"))),
      apply_QCD_LO(string2bool(ctx.get("VJetsReweighting_do_QCD_LO"))),
      apply_QCD_NLO(string2bool(ctx.get("VJetsReweighting_do_QCD_NLO"))),
      apply_QCD_NNLO(string2bool(ctx.get("VJetsReweighting_do_QCD_NNLO"))),
      h_weight_applied(ctx.declare_event_output<float>(weight_name + "_applied")),
      h_weight_EWK(ctx.declare_event_output<float>(weight_name + "_EWK")),
      h_weight_QCD_EWK(ctx.declare_event_output<float>(weight_name + "_QCD_EWK")),
      h_weight_QCD_LO(ctx.declare_event_output<float>(weight_name + "_QCD_LO")),
      h_weight_QCD_NLO(ctx.declare_event_output<float>(weight_name + "_QCD_NLO")),
      h_weight_QCD_NNLO(ctx.declare_event_output<float>(weight_name + "_QCD_NNLO")),
      h_weight_QCD_d1K_NLO(ctx.declare_event_output<float>(weight_name + "_QCD_NLO_d1K")),
      h_weight_QCD_d2K_NLO(ctx.declare_event_output<float>(weight_name + "_QCD_NLO_d2K")),
      h_weight_QCD_d3K_NLO(ctx.declare_event_output<float>(weight_name + "_QCD_NLO_d3K")),
      h_weight_EWK_d1K(ctx.declare_event_output<float>(weight_name + "_EWK_d1K")),
      h_weight_EWK_d2K(ctx.declare_event_output<float>(weight_name + "_EWK_d2K")),
      h_weight_EWK_d3K(ctx.declare_event_output<float>(weight_name + "_EWK_d3K")) {
  // if ((apply_QCD_EWK && (apply_EWK || apply_QCD_NLO)) ||
  //    (apply_QCD_NNLO && !(apply_QCD_EWK || (apply_QCD_NLO)))) {
  //  throw invalid_argument("VJetsReweighting: You are not allowed to use the specified "
  //                        "combination of correction scale factors.");
  // }
  load_correction_files();
}

void VJetsReweighting::load_correction_files() {
  string filesDir = (string)getenv("CMSSW_BASE") + "/src/UHH2/AZH/data/ScaleFactors/VJetsCorrections/";

  // Load EWK correction hiostgrams
  // File naming convention
  // vvj : Z(nunu)
  // eej : Z(ll)
  // evj : W(lv)
  // aj  : Gamma

  // Load Z(nunu)
  TFile *file_vvj = new TFile((filesDir + "vvj.root").c_str());
  for (const string &systematic : {"kappa", "d1kappa", "d2kappa", "d3kappa"})
    load_histo(file_vvj, (string)("vvj_pTV_" + systematic + "_EW"));
  for (const string &pqcd_order : {"NLO", "NNLO"}) {
    for (const string &systematic : {"d1K", "d2K", "d3K"}) {
      load_histo(file_vvj, (string)("vvj_pTV_" + systematic + "_" + pqcd_order));
    }
  }
  load_histo(file_vvj, (string)("vvj_pTV_K_LO"));
  load_histo(file_vvj, (string)("vvj_pTV_K_NLO"));
  load_histo(file_vvj, (string)("vvj_pTV_K_NNLO"));
  file_vvj->Close();

  // Load Z(ll)
  TFile *file_eej = new TFile((filesDir + "eej.root").c_str());
  for (const string &systematic : {"kappa", "d1kappa", "d2kappa", "d3kappa"})
    load_histo(file_eej, (string)("eej_pTV_" + systematic + "_EW"));
  for (const string &pqcd_order : {"NLO", "NNLO"}) {
    for (const string &systematic : {"d1K", "d2K", "d3K"}) {
      load_histo(file_eej, (string)("eej_pTV_" + systematic + "_" + pqcd_order));
    }
  }
  load_histo(file_eej, (string)("eej_pTV_K_LO"));
  load_histo(file_eej, (string)("eej_pTV_K_NLO"));
  load_histo(file_eej, (string)("eej_pTV_K_NNLO"));
  file_eej->Close();
  
  // Load (W(lv))
  TFile *file_evj = new TFile((filesDir + "evj.root").c_str());
  for (const string &systematic : {"kappa", "d1kappa", "d2kappa", "d3kappa"})
    load_histo(file_evj, (string)("evj_pTV_" + systematic + "_EW"));
  for (const string &pqcd_order : {"NLO", "NNLO"}) {
    for (const string &systematic : {"d1K", "d2K", "d3K"}) {
      load_histo(file_evj, (string)("evj_pTV_" + systematic + "_" + pqcd_order));
    }
  }
  load_histo(file_evj, (string)("evj_pTV_K_LO"));
  load_histo(file_evj, (string)("evj_pTV_K_NLO"));
  load_histo(file_evj, (string)("evj_pTV_K_NNLO"));
  file_evj->Close();
}

void VJetsReweighting::load_histo(TFile *file, const string &name) {
  histos[name].reset((TH1F *)file->Get(name.c_str()));
  histos[name]->SetDirectory(0);
}

double VJetsReweighting::get_gen_v_pt(const Event &event) {
  if (!(is_DYJets || is_WJets  || is_ZJets))
    throw runtime_error("VJetsReweighting::get_gen_v_pt(): Calling this function "
                        "on non-WJets/DYJets sample makes no sense.");
  double pt(-1.);
  bool v_found(false);
  for (const GenParticle &gp : *event.genparticles) {
    if (is_WJets && gp.status() == 22 && abs(gp.pdgId()) == 24) {
      pt = gp.v4().Pt();
      v_found = true;
    } else if (is_DYJets && gp.status() == 22 && abs(gp.pdgId()) == 23) {
      pt = gp.v4().Pt();
      v_found = true;
    } else if (is_ZJets && gp.status() == 22 && abs(gp.pdgId()) == 23) {
      pt = gp.v4().Pt();
      v_found = true;
    }
  }
  if (!v_found) {
    int n_status23_leptons(0);
    GenParticle d1, d2; // daughters of V boson
    for (const GenParticle &gp : *event.genparticles) {
      if (gp.status() == 23 && abs(gp.pdgId()) >= 11 && abs(gp.pdgId()) <= 16) {
        n_status23_leptons++;
        if (gp.pdgId() > 0)
          d1 = gp;
        else
          d2 = gp;
      }
    }
    if (n_status23_leptons != 2)
      throw runtime_error("VJetsReweighting::get_gen_v_pt(): Did not find exactly "
                          "two V daughter candidates.");
    pt = (d1.v4() + d2.v4()).Pt();
  }
  return pt;
}

double VJetsReweighting::evaluate(const string &name, const double pt) {
  const int firstBin = 1;
  const int lastBin = histos[name]->GetNbinsX();
  const double h_min = histos[name]->GetBinCenter(firstBin) - 0.5 * histos[name]->GetBinWidth(firstBin);
  const double h_max = histos[name]->GetBinCenter(lastBin) + 0.5 * histos[name]->GetBinWidth(lastBin);
  double pt_for_eval = pt;
  pt_for_eval = (pt_for_eval > h_min) ? pt_for_eval : h_min + 0.001;
  pt_for_eval = (pt_for_eval < h_max) ? pt_for_eval : h_max - 0.001;
  return histos[name]->GetBinContent(histos[name]->FindBin(pt_for_eval));
}

bool VJetsReweighting::process(Event &event) {
  float weight_applied(1.);
  float weight_EWK(1.);
  float weight_QCD_EWK(1.);
  float weight_QCD_LO(1.);
  float weight_QCD_NLO(1.);
  float weight_QCD_NNLO(1.);
  float weight_QCD_d1K_NLO(1.);
  float weight_QCD_d2K_NLO(1.);
  float weight_QCD_d3K_NLO(1.);
  float weight_EWK_d1K(1.);
  float weight_EWK_d2K(1.);
  float weight_EWK_d3K(1.);

  if (!event.isRealData && (is_WJets || is_DYJets || is_ZJets)) {
    double pt = get_gen_v_pt(event);
    string file("");
    if (is_WJets)
      file = "evj";
    if (is_DYJets)
      file = "eej";
    if (is_ZJets)
      file = "vvj";

    // QCD EWK
    if (apply_QCD_EWK) {
      throw logic_error("QCD EWK not implemented yet.");
      // weight_QCD_EWK = evaluate(file + "", pt);
      // weight_applied *= weight_QCD_EWK;
    }
    // EWK
    if (apply_EWK) {
      weight_EWK = evaluate(file + "_pTV_kappa_EW", pt);
      weight_EWK_d1K = evaluate(file + "_pTV_d1kappa_EW", pt);
      weight_EWK_d2K = evaluate(file + "_pTV_d2kappa_EW", pt);
      weight_EWK_d3K = evaluate(file + "_pTV_d3kappa_EW", pt);
      weight_applied *= 1 + weight_EWK;
    }
    // QCD LO
    if (apply_QCD_LO) {
      weight_QCD_LO = evaluate(file + "_pTV_K_LO", pt);
      weight_applied *= weight_QCD_LO;
    }
    // QCD NLO
    if (apply_QCD_NLO) {
      weight_QCD_NLO = evaluate(file + "_pTV_K_NLO", pt);
      weight_QCD_d1K_NLO = evaluate(file + "_pTV_d1K_NLO", pt);
      weight_QCD_d2K_NLO = evaluate(file + "_pTV_d2K_NLO", pt);
      weight_QCD_d3K_NLO = evaluate(file + "_pTV_d3K_NLO", pt);
      weight_applied *= weight_QCD_NLO;
    }
    // QCD NNLO
    if (apply_QCD_NNLO) {
      // throw logic_error("QCD NNLO not implemented yet.");
      weight_QCD_NNLO = evaluate(file + "_pTV_K_NNLO", pt);
      weight_applied *= weight_QCD_NNLO;
    }
  }

  event.weight *= weight_applied;
  event.set(h_weight_applied, weight_applied);
  event.set(h_weight_EWK, weight_EWK);
  event.set(h_weight_QCD_EWK, weight_QCD_EWK);
  event.set(h_weight_QCD_LO, weight_QCD_LO);
  event.set(h_weight_QCD_NLO, weight_QCD_NLO);
  event.set(h_weight_QCD_NNLO, weight_QCD_NNLO);
  // Systematic uncertainties
  event.set(h_weight_QCD_d1K_NLO, weight_QCD_d1K_NLO);
  event.set(h_weight_QCD_d2K_NLO, weight_QCD_d2K_NLO);
  event.set(h_weight_QCD_d3K_NLO, weight_QCD_d3K_NLO);
  event.set(h_weight_EWK_d1K, weight_EWK_d1K);
  event.set(h_weight_EWK_d2K, weight_EWK_d2K);
  event.set(h_weight_EWK_d3K, weight_EWK_d3K);

  return true;
}