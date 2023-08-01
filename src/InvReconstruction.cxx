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
#include "UHH2/AZH/include/HiggsReco.h"
#include "UHH2/AZH/include/AtoZHHists.h"
#include "UHH2/AZH/include/Utils.h"

using namespace std;
using namespace uhh2;

class InvReconstruction : public AnalysisModule{
  public:
  explicit InvReconstruction(Context & ctx);
  virtual bool process(Event& event) override;


  private:

  //Methods
  bool AssignRegion(Event& event);
  bool HasGenLepton(Event& event);
  double DeltaPhi(Event& event);

  //Handles
  Event::Handle<int> handle_region;
  Event::Handle<double> handle_sum_pt;
  Event::Handle<double> handle_delta;
  Event::Handle<double> handle_weight;
  Event::Handle<double> handle_H_m;
  Event::Handle<double> handle_met;
  Event::Handle<double> handle_A_mt;
  Event::Handle<double> handle_H_mt;
  Event::Handle<double> handle_mt_diff;
  Event::Handle<double> handle_m_diff;
    
  //Histograms
  unique_ptr<Hists> h_base;
  unique_ptr<Hists> h_met;
  unique_ptr<Hists> h_weight;
  unique_ptr<Hists> h_delta;
  unique_ptr<Hists> h_btag;
  unique_ptr<Hists> h_veto;
  unique_ptr<Hists> h_6jets;
  unique_ptr<Hists> h_fake;
  unique_ptr<Hists> h_missed;    

  //Scale Factors  
  std::unique_ptr<AnalysisModule> sf_btagging;
  //std::unique_ptr<AnalysisModule> sf_leptons;
    
  //Selection Modules
  std::unique_ptr<Selection> s_njet_six;
  std::unique_ptr<Selection> s_bjet_one;
  std::unique_ptr<Selection> s_bjet_two;

  //Other Fields
  std::unique_ptr<HiggsReconstructor> higgs_reconstructor;
  const JetId bmedium = BTag(BTag::DEEPJET, BTag::WP_MEDIUM);

};


InvReconstruction::InvReconstruction(Context& ctx){

  sf_btagging.reset(new MCBTagScaleFactor(ctx, BTag::DEEPJET, BTag::WP_MEDIUM, "jets", "mujets", "incl","BTagMCEffFile"));
  //sf_leptons.reset(new LeptonScaleFactors(ctx));
      
  s_njet_six.reset(new NJetSelection(6));
  s_bjet_one.reset(new NJetSelection(1,1,bmedium));
  s_bjet_two.reset(new NJetSelection(2,-1,bmedium));

  h_base.reset(new SimpleHist(ctx, "CutFlow_Baseline"));
  h_met.reset(new SimpleHist(ctx, "CutFlow_MET>170"));
  h_weight.reset(new SimpleHist(ctx, "CutFlow_MET>170*"));
  h_delta.reset(new SimpleHist(ctx, "CutFlow_deltaphi"));
  h_btag.reset(new SimpleHist(ctx, "CutFlow_TwoB"));
  h_veto.reset(new SimpleHist(ctx, "CutFlow_LeptonVeto"));
  h_6jets.reset(new SimpleHist(ctx, "CutFlow_6Jets"));
  h_fake.reset(new SimpleHist(ctx, "faked_leptons"));
  h_missed.reset(new SimpleHist(ctx, "missed_leptons"));

  higgs_reconstructor.reset(new HiggsReconstructor(ctx));

  handle_weight = ctx.get_handle<double>("event_weight");
  handle_region = ctx.declare_event_output<int>("region");
  handle_sum_pt = ctx.declare_event_output<double>("HT");
  handle_delta = ctx.declare_event_output<double>("delta_phi");
  handle_met = ctx.declare_event_output<double>("MET");
  handle_A_mt = ctx.declare_event_output<double>("mt_A");
  handle_H_mt = ctx.declare_event_output<double>("mt_H");
  handle_H_m = ctx.declare_event_output<double>("m_H");
  handle_mt_diff = ctx.declare_event_output<double>("delta_mt");
  handle_m_diff = ctx.declare_event_output<double>("delta_m");

}


bool InvReconstruction::process(Event& event){

  event.weight = event.get(handle_weight);

  //BTag scale factors
  sf_btagging->process(event);

  //Sanity check for lepton sf
  int leptons = (*event.electrons).size() + (*event.muons).size();
  if(leptons==0&&HasGenLepton(event)) h_missed->fill(event);
  if(leptons!=0&&(!HasGenLepton(event))) h_fake->fill(event);

  //Global cuts
  h_base->fill(event);
  if (event.met->pt() < 170) return false;
  h_met->fill(event);
  if (event.weight > 10) return false;
  h_weight->fill(event);
  if (DeltaPhi(event) < 0.5) return false;
  h_delta->fill(event);

  //Rest of cutflow
  int lep = (*event.electrons).size() + (*event.muons).size();
  if (lep==0) h_veto->fill(event);
  if (lep==0 && s_bjet_two->passes(event)) h_btag->fill(event);
  if (lep==0 && s_bjet_two->passes(event) && s_njet_six->passes(event)) h_6jets->fill(event);
   
  //Assign regions
  bool valid_region = AssignRegion(event);
  if(!valid_region) return false;

  //Reconstruct A and H
  higgs_reconstructor->process(event);
  event.set(handle_A_mt,higgs_reconstructor->GetTransMassA());
  event.set(handle_H_mt,higgs_reconstructor->GetTransMassH());
  event.set(handle_H_m,higgs_reconstructor->GetMassH());
  event.set(handle_mt_diff,higgs_reconstructor->GetTransMassDiff());
  event.set(handle_m_diff,higgs_reconstructor->GetMassDiff());  

  //Setting handles
  double HT = 0;
  for(Jet jet: *event.jets){HT += jet.pt();}

  event.set(handle_sum_pt,HT);
  event.set(handle_delta,DeltaPhi(event));
  event.set(handle_met,event.met->pt());
  event.set(handle_weight, event.weight);

  return true;
}

bool InvReconstruction::AssignRegion(Event& event){
  int leptons = (*event.electrons).size() + (*event.muons).size();
  int btags = 0;
  if (s_bjet_one->passes(event)) btags = 1;
  if (s_bjet_two->passes(event)) btags = 2;
  
  int jets = 5;
  if (s_njet_six->passes(event)) jets = 6;

  array<int,3> lbj = {leptons,btags,jets};
  Region region = Region::Invalid;

  if(lbj==array<int,3>{0,2,6}) region = Region::SR_6J;    //Signal region
  if(lbj==array<int,3>{0,2,5}) region = Region::SR_5J; 
  if(lbj==array<int,3>{0,1,6}) region = Region::IR_1B_6J; //Invisible control region
  if(lbj==array<int,3>{0,1,5}) region = Region::IR_1B_5J;
  if(lbj==array<int,3>{0,0,6}) region = Region::IR_0B_6J;
  if(lbj==array<int,3>{0,0,5}) region = Region::IR_0B_5J;
  if(lbj==array<int,3>{1,2,6}) region = Region::LR_2B_6J; //Leptonic control region
  if(lbj==array<int,3>{1,2,5}) region = Region::LR_2B_5J;
  if(lbj==array<int,3>{1,1,6}) region = Region::LR_1B_6J;
  if(lbj==array<int,3>{1,1,5}) region = Region::LR_1B_5J;
  if(lbj==array<int,3>{1,0,6}) region = Region::LR_0B_6J;
  if(lbj==array<int,3>{1,0,5}) region = Region::LR_0B_5J;

  if(region==Region::Invalid) return false;

  event.set(handle_region, (int) region);
  return true;
}

bool InvReconstruction::HasGenLepton(Event& event){
  for(GenParticle P : *event.genparticles){
      if (abs(P.pdgId()) == 11) return true;
      if (abs(P.pdgId()) == 13) return true;
  }
  return false;
}

double InvReconstruction::DeltaPhi(Event& event){
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

UHH2_REGISTER_ANALYSIS_MODULE(InvReconstruction)