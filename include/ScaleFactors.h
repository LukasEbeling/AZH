#pragma once
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Electron.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Muon.h"

#include "UHH2/common/include/YearRunSwitchers.h"

#include <TH2F.h>

class ScaleFactors2016preVFP : public uhh2::AnalysisModule {
public:
  explicit ScaleFactors2016preVFP(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

private:
  std::unique_ptr<AnalysisModule> mu_sf_reco,mu_sf_id, mu_sf_iso, ele_sf_id, ele_sf_reco;
};

class ScaleFactors2016postVFP : public uhh2::AnalysisModule {
public:
  explicit ScaleFactors2016postVFP(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

private:
  std::unique_ptr<AnalysisModule> mu_sf_reco,mu_sf_id, mu_sf_iso, ele_sf_id, ele_sf_reco;
};

class ScaleFactors2017 : public uhh2::AnalysisModule {
public:
  explicit ScaleFactors2017(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

private:
  std::unique_ptr<AnalysisModule> mu_sf_reco,mu_sf_id, mu_sf_iso, ele_sf_id, ele_sf_reco;
};

class ScaleFactors2018 : public uhh2::AnalysisModule {
public:
  explicit ScaleFactors2018(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

private:
  std::unique_ptr<AnalysisModule> mu_sf_reco,mu_sf_id, mu_sf_iso, ele_sf_id, ele_sf_reco;
};

class LeptonScaleFactors : public uhh2::AnalysisModule {
public:
  explicit LeptonScaleFactors(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event);

private:
  std::unique_ptr<YearSwitcher> m_sf_lepton;
};

class PUIDScaleFactors : public uhh2::AnalysisModule {
public:
  explicit PUIDScaleFactors(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event);
  double load_sf(const Jet &jet, std::string variation);

private:
  Year year;
  bool is_mc;
  uhh2::Event::Handle<float> handle_sf_puid, handle_sf_puid_up, handle_sf_puid_down;
  std::unique_ptr<TFile> file_sf;
};

class L1PrefiringWeight : public uhh2::AnalysisModule {
public:
  explicit L1PrefiringWeight(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event);

private:
  Year year;
  int syst_direction;
};

class DiEleTriggerScaleFactors : public uhh2::AnalysisModule {
public:
  explicit DiEleTriggerScaleFactors(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

private:
  Year year;
  int syst_direction;
  uhh2::Event::Handle<float> h_ee_trigger_sf, h_ee_trigger_sf_up, h_ee_trigger_sf_dn;
  uhh2::Event::Handle<std::vector<Electron>> handle_tight_electrons;
};

class DiMuTriggerScaleFactors : public uhh2::AnalysisModule {
public:
  explicit DiMuTriggerScaleFactors(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

private:
  Year year;
  int syst_direction;
  uhh2::Event::Handle<float> h_mumu_trigger_sf, h_mumu_trigger_sf_up, h_mumu_trigger_sf_dn;
  uhh2::Event::Handle<std::vector<Muon>> handle_tight_muons;
};

class EleMuTriggerScaleFactors : public uhh2::AnalysisModule {
public:
  explicit EleMuTriggerScaleFactors(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

private:
  Year year;
  int syst_direction;
  uhh2::Event::Handle<float> h_emu_trigger_sf, h_emu_trigger_sf_up, h_emu_trigger_sf_dn;
};

class InitTriggerScaleFactors : public uhh2::AnalysisModule {
public:
  explicit InitTriggerScaleFactors(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event);

private:
  uhh2::Event::Handle<int> handle_channel;
  std::vector<uhh2::Event::Handle<float>> m_handles;
};

class TriggerScaleFactors : public uhh2::AnalysisModule {
public:
  explicit TriggerScaleFactors(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event);

private:
  uhh2::Event::Handle<int> handle_channel;
  std::unique_ptr<uhh2::AnalysisModule> m_sf_ee_trigger;
  std::unique_ptr<uhh2::AnalysisModule> m_sf_mumu_trigger;
  std::unique_ptr<uhh2::AnalysisModule> m_sf_emu_trigger;
  std::unique_ptr<uhh2::AnalysisModule> m_sf_dummy;
};

//____________________________________________________________________________________________________
// Copy of https://github.com/MatthiesC/HighPtSingleTop/blob/master/include/TheoryCorrections.h
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting

class TopPtReweighting : public uhh2::AnalysisModule {
public:
  TopPtReweighting(uhh2::Context &ctx, const bool apply = true);
  virtual bool process(uhh2::Event &event) override;

private:
  uhh2::Event::Handle<float> h_weight_nominal;
  uhh2::Event::Handle<float> h_weight_a_up;
  uhh2::Event::Handle<float> h_weight_a_down;
  uhh2::Event::Handle<float> h_weight_b_up;
  uhh2::Event::Handle<float> h_weight_b_down;
  uhh2::Event::Handle<float> h_weight_applied;
  bool proc_is_TT;
  void set_dummy_weights(uhh2::Event &event);
  const float fDummyWeight = 1.0f;
  const bool fApply;
  enum class TopPtVariation {
    nominal,
    a_up,
    a_down,
    b_up,
    b_down,
  };
  TopPtVariation applied_variation;
  const bool fPtCutOff_b = false;
  const float fPtCutOff = 500.;
  const float fA = 0.0615;
  const float fA_up = fA * 1.5;
  const float fA_down = fA * 0.5;
  const float fB = -0.0005;
  const float fB_up = fB * 1.5;
  const float fB_down = fB * 0.5;
};

/*
Original version taken from https://github.com/UHH2/VHResonances/blob/master/Analysis/python/PlotNLOCorrections.py
Example module can be found here: https://github.com/UHH2/VHResonances/blob/master/src/HiggsToWWModules.cxx
*/
class VJetsReweighting : public uhh2::AnalysisModule {
public:
  explicit VJetsReweighting(uhh2::Context &ctx, const std::string &weight_name = "weight_vjets");
  virtual bool process(uhh2::Event &event) override;

private:
  std::unordered_map<std::string, std::unique_ptr<TH1F>> histos;

  void load_correction_files();
  void load_histo(TFile *file, const std::string &name);
  double get_gen_v_pt(const uhh2::Event &event);
  double evaluate(const std::string &name, const double pt);

  const bool is_WJets;
  const bool is_DYJets;
  const bool apply_EWK;
  const bool apply_QCD_EWK;
  const bool apply_QCD_LO;
  const bool apply_QCD_NLO;
  const bool apply_QCD_NNLO;

  const uhh2::Event::Handle<float> h_weight_applied;
  const uhh2::Event::Handle<float> h_weight_EWK;
  const uhh2::Event::Handle<float> h_weight_QCD_EWK;
  const uhh2::Event::Handle<float> h_weight_QCD_LO;
  const uhh2::Event::Handle<float> h_weight_QCD_NLO;
  const uhh2::Event::Handle<float> h_weight_QCD_NNLO;
  // Systematic uncertainties
  const uhh2::Event::Handle<float> h_weight_QCD_d1K_NLO;
  const uhh2::Event::Handle<float> h_weight_QCD_d2K_NLO;
  const uhh2::Event::Handle<float> h_weight_QCD_d3K_NLO;
  const uhh2::Event::Handle<float> h_weight_EWK_d1K;
  const uhh2::Event::Handle<float> h_weight_EWK_d2K;
  const uhh2::Event::Handle<float> h_weight_EWK_d3K;
};