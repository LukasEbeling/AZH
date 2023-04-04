#pragma once
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/YearRunSwitchers.h"

#include <TH2F.h>

class ScaleFactors2016preVFP: public uhh2::AnalysisModule {
 public:
  explicit ScaleFactors2016preVFP(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;

 private:
  std::unique_ptr<AnalysisModule> mu_sf_id, mu_sf_iso, ele_sf_id, ele_sf_reco;
};


class ScaleFactors2016postVFP: public uhh2::AnalysisModule {
 public:
  explicit ScaleFactors2016postVFP(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;

 private:
  std::unique_ptr<AnalysisModule> mu_sf_id, mu_sf_iso, ele_sf_id, ele_sf_reco;
};


class ScaleFactors2017: public uhh2::AnalysisModule {
 public:
  explicit ScaleFactors2017(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;

 private:
  std::unique_ptr<AnalysisModule> mu_sf_id, mu_sf_iso, ele_sf_id, ele_sf_reco;
};


class ScaleFactors2018: public uhh2::AnalysisModule {
 public:
  explicit ScaleFactors2018(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;

 private:
  std::unique_ptr<AnalysisModule> mu_sf_id, mu_sf_iso, ele_sf_id, ele_sf_reco;
};


class LeptonScaleFactors: public uhh2::AnalysisModule {
 public:
  explicit LeptonScaleFactors(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event);

 private:
  std::unique_ptr<YearSwitcher> m_sf_lepton;
};


class L1PrefiringWeight: public uhh2::AnalysisModule {
  public:
    explicit L1PrefiringWeight(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & event);

  private:
    Year year;
    int syst_direction;
};


class DiEleTriggerScaleFactors: public uhh2::AnalysisModule {
  public:
    explicit DiEleTriggerScaleFactors(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & event) override;

  private:
    Year year;
    int syst_direction;
    uhh2::Event::Handle<float> h_ee_trigger_sf, h_ee_trigger_sf_up, h_ee_trigger_sf_dn;
};


class DiMuTriggerScaleFactors: public uhh2::AnalysisModule {
  public:
    explicit DiMuTriggerScaleFactors(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & event) override;

  private:
    Year year;
    int syst_direction;
    uhh2::Event::Handle<float> h_mumu_trigger_sf, h_mumu_trigger_sf_up, h_mumu_trigger_sf_dn;
};


class EleMuTriggerScaleFactors: public uhh2::AnalysisModule {
  public:
    explicit EleMuTriggerScaleFactors(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & event) override;

  private:
    Year year;
    int syst_direction;
    uhh2::Event::Handle<float> h_emu_trigger_sf, h_emu_trigger_sf_up, h_emu_trigger_sf_dn;
};


class InitTriggerScaleFactors: public uhh2::AnalysisModule {
  public:
    explicit InitTriggerScaleFactors(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & event);

  private:
    uhh2::Event::Handle<int> handle_channel;
    std::vector<uhh2::Event::Handle<float>> m_handles;
};


class TriggerScaleFactors: public uhh2::AnalysisModule {
  public:
    explicit TriggerScaleFactors(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & event);

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

class TopPtReweighting: public uhh2::AnalysisModule {
public:
  TopPtReweighting(uhh2::Context & ctx, const bool apply = true);
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<float> h_weight_nominal;
  uhh2::Event::Handle<float> h_weight_a_up;
  uhh2::Event::Handle<float> h_weight_a_down;
  uhh2::Event::Handle<float> h_weight_b_up;
  uhh2::Event::Handle<float> h_weight_b_down;
  uhh2::Event::Handle<float> h_weight_applied;
  bool proc_is_TT;
  void set_dummy_weights(uhh2::Event & event);
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
  const float fA_up = fA*1.5;
  const float fA_down = fA*0.5;
  const float fB = -0.0005;
  const float fB_up = fB*1.5;
  const float fB_down = fB*0.5;
};


/*
Taken from https://github.com/UHH2/VHResonances/blob/master/Analysis/python/PlotNLOCorrections.py
- EWK corrections are available for W+jets, Z+jets, gamma+jets samples in the "merged_kfactors_*.root" files under the name "kfactor_monojet_ewk"
- QCD NLO corrections are available for W+jets, Z+jets, gamma+jets in the same files under the name "kfactor_monojet_qcd". Those are calculated for 2016 samples.
- 2017 version of QCD NLO corrections are available for Z+jets (ll + nunu cases) in the "kfac_*_filter" files.
- QCD NNLO corrections are in the "lindert_qcd_nnlo_sf" file with the following convention:
    - eej -> Z(ll) +jets
    - vvj -> Z(nunu) +jets
    - evj -> W +jets
    - aj -> gamma +jets
- QCD NNLO corrections need to be applied on top of EWK corrections for NLO samples and on top of EWK + QCD NLO corrections for LO samples.
- According to Andreas "I do not apply the NNLO corrections. I have not seen any evidence that they actually improve data/MC agreement. I do not trust them."
- For W+Jets @LO for 2017 and 2018: wjet_dress_monojet or wjet_dress_inclusive in "2017_gen_v_pt_qcd_sf.root"
- In the "merged_kfactors_*.root" file, for Z and W + jets, the qcd_ewk histograms are also present: qcd_ewk = QCD * EWK
- taken from https://github.com/bu-cms/bucoffea/tree/master/bucoffea/data/sf/theory
- relative to those studies https://arxiv.org/abs/1705.04664
Example module can be found here: https://github.com/UHH2/VHResonances/blob/master/src/HiggsToWWModules.cxx
*/

class VJetsReweighting: public uhh2::AnalysisModule {
  public:
   explicit VJetsReweighting(uhh2::Context & ctx, const std::string& weight_name="weight_vjets");
   virtual bool process(uhh2::Event & event) override;

  private:
   std::unordered_map<std::string, std::unique_ptr<TH1F>> histos;

   void load_histo(TFile* file, const std::string& name, const std::string& histName);
   double get_v_pt(const uhh2::Event & event);
   double evaluate(const std::string& name, const double pt);

   const bool is_2016_nonUL;
   const bool is_WJets;
   const bool is_DYJets;
   const bool apply_EWK;
   const bool apply_QCD_EWK;
   const bool apply_QCD_NLO;
   const bool apply_QCD_NNLO;

   const uhh2::Event::Handle<float> h_weight_applied;
   const uhh2::Event::Handle<float> h_weight_EWK;
   const uhh2::Event::Handle<float> h_weight_QCD_EWK;
   const uhh2::Event::Handle<float> h_weight_QCD_NLO;
   const uhh2::Event::Handle<float> h_weight_QCD_NNLO;
};

