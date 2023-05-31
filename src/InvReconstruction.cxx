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

  //Scale Factors  
  std::unique_ptr<AnalysisModule> sf_btagging;
  
  //Selection Modules
  std::unique_ptr<Selection> s_bjet_two;
};

InvReconstruction::InvReconstruction(Context& ctx){
  handle_region = ctx.declare_event_output<int>("region");
  handle_sum_pt = ctx.declare_event_output<double>("HT");
  handle_weight = ctx.get_handle<double>("event_weight");
  handle_delta = ctx.declare_event_output<double>("delta_phi");
  handle_leptons = ctx.declare_event_output<int>("num_leptons");

  higgs_reconstructor.reset(new HiggsReconstructor(ctx));

  run_btag_sf = ctx.has("BTagMCEffFile");
  if (run_btag_sf){
    sf_btagging.reset(new MCBTagScaleFactor(ctx, BTag::DEEPJET, BTag::WP_MEDIUM, "jets", "mujets", "incl","BTagMCEffFile"));
  }

  s_bjet_two.reset(new NJetSelection(2,-1,bmedium));
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

  int leptons = (*event.electrons).size() + (*event.muons).size();
  event.set(handle_leptons,leptons);
  
  Region region = Region::Invalid;
  if (met_high && delta_high && two_b && leptons == 0) region = Region::SignalRegion;
  else if (met_high && delta_high && two_b && leptons == 1) region = Region::CRTTBar;
  else if (met_high && !delta_high && two_b && leptons == 0) region = Region::CRQCD;
  else if (met_high && delta_high && true && leptons == 2) region = Region::CRZJets;
  else if (!met_high && delta_high && two_b && leptons == 0) region = Region::CRMet;
  else if (true && delta_high && two_b && leptons == 1) region = Region::CRTTBar2;
  else return false;

  event.set(handle_region, (int) region);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(InvReconstruction)