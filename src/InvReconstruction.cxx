#include <iostream>
#include <memory>
#include <cmath>

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/PSWeights.h"

#include "UHH2/AZH/include/NormalisationTools.h"
#include "UHH2/AZH/include/ScaleFactors.h"
#include "UHH2/AZH/include/Utils.h"//for region
#include "UHH2/AZH/include/HiggsReco.h"
#include "UHH2/AZH/include/AtoZHHists.h"
#include "UHH2/AZH/include/METTriggers.h"


using namespace std;
using namespace uhh2;

class InvReconstruction : public AnalysisModule{
  public:
  explicit InvReconstruction(Context& ctx);
  virtual bool process(Event& event) override;

  private:
  unique_ptr<HiggsReconstructor> higgs_reconstructor;
  bool run_btag_sf;
  const JetId bmedium = BTag(BTag::DEEPJET, BTag::WP_MEDIUM);
  
  //Methods
  bool assign_region(Event& event);
  double delta_phi(Event& event);

  //Handles
  Event::Handle<int> handle_region;
  Event::Handle<double> handle_sum_pt;
  Event::Handle<double> handle_weight;
  Event::Handle<double> handle_delta;
  Event::Handle<int> handle_leptons;
  Event::Handle<int> handle_mtrigger;

  //Histograms
  unique_ptr<Hists> h_base;
  unique_ptr<Hists> h_met;
  unique_ptr<Hists> h_delta;
  unique_ptr<Hists> h_btag;
  unique_ptr<Hists> h_veto;
  unique_ptr<Hists> h_weight;

  //Scale Factors  
  std::unique_ptr<AnalysisModule> sf_btagging;
  
  //Selection Modules
  std::unique_ptr<Selection> s_bjet_two;
  std::unique_ptr<Selection> s_bjet_none;
  std::unique_ptr<METTriggers> s_met_trigger;
};

InvReconstruction::InvReconstruction(Context& ctx){
  handle_region = ctx.declare_event_output<int>("region");
  handle_sum_pt = ctx.declare_event_output<double>("HT");
  handle_weight = ctx.get_handle<double>("event_weight");
  handle_delta = ctx.declare_event_output<double>("delta_phi");
  handle_leptons = ctx.declare_event_output<int>("num_leptons");
  handle_mtrigger = ctx.declare_event_output<int>("met_triggered");

  h_base.reset(new RecoHists(ctx, "CutFlow_Baseline"));
  h_met.reset(new RecoHists(ctx, "CutFlow_MET>170"));
  h_delta.reset(new RecoHists(ctx, "CutFlow_deltaphi"));
  h_btag.reset(new RecoHists(ctx, "CutFlow_TwoB"));
  h_veto.reset(new RecoHists(ctx, "CutFlow_LeptonVeto"));
  h_weight.reset(new RecoHists(ctx, "CutFlow_Weights"));


  higgs_reconstructor.reset(new HiggsReconstructor(ctx));

  run_btag_sf = ctx.has("BTagMCEffFile");
  if (run_btag_sf){
    sf_btagging.reset(new MCBTagScaleFactor(ctx, BTag::DEEPJET, BTag::WP_MEDIUM, "jets", "mujets", "incl","BTagMCEffFile"));
  }

  s_bjet_two.reset(new NJetSelection(2,-1,bmedium));
  s_bjet_none.reset(new NJetSelection(-1,0,bmedium));
  s_met_trigger.reset(new METTriggers(ctx));
}

bool InvReconstruction::process(Event& event){
  event.weight = event.get(handle_weight);

  if (run_btag_sf) sf_btagging->process(event);

  bool valid_region = assign_region(event);
  if(!valid_region) return false;

  event.set(handle_delta,delta_phi(event));
  higgs_reconstructor->process(event);

  double HT = 0;
  for(Jet jet: *event.jets){HT += jet.pt();}
  event.set(handle_sum_pt,HT);

  int triggered = 0;
  if (s_met_trigger->passes(event)) triggered = 1;
  event.set(handle_mtrigger,triggered);

  event.set(handle_weight, event.weight);

  return true;
}

double InvReconstruction::delta_phi(Event& event){
  double phi_m = event.met->phi();
  double min_diff = M_PI;

  for (Jet jet: *event.jets){
    if(jet.pt()<30) continue;
    double phi_j = jet.phi();
    double delta = abs(phi_m - phi_j);
    if(delta > M_PI){delta = 2*M_PI - delta;}
    if(delta < min_diff) min_diff = delta;
  }

  return min_diff;
}

bool InvReconstruction::assign_region(Event& event){
  bool met_high = event.met->pt() > 170;
  bool delta_high = delta_phi(event) > 0.5;
  bool two_b = s_bjet_two->passes(event);
  bool no_b = s_bjet_none->passes(event);

  int leptons = (*event.electrons).size() + (*event.muons).size();
  event.set(handle_leptons,leptons);
  
  //h_base -> fill(event);
  //if (met_high) h_met -> fill(event);
  //if (met_high && delta_high) h_delta -> fill(event);
  //if (met_high && delta_high && two_b) h_btag -> fill(event);
  //if (met_high && delta_high && two_b && leptons == 0) h_veto -> fill(event);
  //if (met_high && delta_high && two_b && leptons == 0 && event.weight <= 10) h_weight -> fill(event);

  h_base -> fill(event);
  if (met_high) h_met -> fill(event);
  if (met_high && event.weight <= 10) h_weight -> fill(event);
  if (met_high && event.weight <= 10 && delta_high) h_delta -> fill(event);
  if (met_high && event.weight <= 10 && delta_high && two_b) h_btag -> fill(event);
  if (met_high && event.weight <= 10 && delta_high && two_b && leptons == 0) h_veto -> fill(event);


  Region region = Region::Invalid;
  if (met_high && delta_high && two_b && leptons == 0 && event.weight < 10) region = Region::SR;
  else if (met_high && delta_high && two_b && leptons == 1 && event.weight < 10) region = Region::CR_1L;
  //else if (met_high && !delta_high && two_b && leptons == 0) region = Region::CR_lowdelta;
  else if (met_high && delta_high && no_b && leptons == 0 && event.weight < 10) region = Region::CR_0B;
  //else if (met_high && delta_high && no_b && leptons == 2) region = Region::CR_0B_2L;
  else if (!met_high && delta_high && two_b && leptons == 0) region = Region::CR_lowmet;
  //else if (true && delta_high && two_b && leptons == 1) region = Region::CR_1L_anymet;
  else return false;

  event.set(handle_region, (int) region);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(InvReconstruction)