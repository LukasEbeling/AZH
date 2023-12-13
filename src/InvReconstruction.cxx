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

  bool is_mc = false;

  //Methods
  bool AssignRegion(Event& event);
  bool HasGenLepton(Event& event);
  double DeltaPhi(Event& event);
  double DeltaEta(Event& event);
  double get_random(Event& event);

  //Handles
  Event::Handle<int> handle_region;
  Event::Handle<int> handle_node;
  Event::Handle<double> handle_HT;
  Event::Handle<double> handle_delta_phi;
  Event::Handle<double> handle_delta_eta;
  Event::Handle<double> handle_weight;
  Event::Handle<double> handle_H_m;
  Event::Handle<double> handle_met;
  Event::Handle<double> handle_A_mt;
  Event::Handle<double> handle_H_mt;
  Event::Handle<double> handle_mt_diff;
  Event::Handle<double> handle_m_diff;
  Event::Handle<double> handle_score;
  Event::Handle<int> handle_num_l;
  Event::Handle<int> handle_num_b;
  Event::Handle<int> handle_num_j;
  Event::Handle<int> handle_kfold;
    
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
  std::unique_ptr<Selection> s_njet_sixplus;
  std::unique_ptr<Selection> s_njet_sevenplus;
  std::unique_ptr<Selection> s_bjet_none;
  std::unique_ptr<Selection> s_bjet_one;
  std::unique_ptr<Selection> s_bjet_two;

  //Other Fields
  std::unique_ptr<HiggsReconstructor> higgs_reconstructor;
  const JetId bmedium = BTag(BTag::DEEPJET, BTag::WP_MEDIUM);

};


InvReconstruction::InvReconstruction(Context& ctx){
  is_mc = ctx.get("dataset_type") == "MC";

  sf_btagging.reset(new MCBTagScaleFactor(ctx, BTag::DEEPJET, BTag::WP_MEDIUM, "jets", "mujets", "incl","BTagMCEffFile"));
  //sf_leptons.reset(new LeptonScaleFactors(ctx));
  
  s_njet_six.reset(new NJetSelection(6,6));
  s_njet_sixplus.reset(new NJetSelection(6));
  s_njet_sevenplus.reset(new NJetSelection(7));
  s_bjet_none.reset(new NJetSelection(0,-1,bmedium));
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
  handle_HT = ctx.declare_event_output<double>("HT");
  handle_met = ctx.declare_event_output<double>("MET");
  handle_A_mt = ctx.declare_event_output<double>("mt_A");
  handle_H_mt = ctx.declare_event_output<double>("mt_H");
  handle_H_m = ctx.declare_event_output<double>("m_H");
  handle_mt_diff = ctx.declare_event_output<double>("delta_mt");
  handle_m_diff = ctx.declare_event_output<double>("delta_m");
  handle_delta_phi = ctx.declare_event_output<double>("delta_phi"); //min delta phi
  handle_delta_eta = ctx.declare_event_output<double>("delta_eta"); //leading two jets
  handle_score = ctx.declare_event_output<double>("score"); //dnn score
  handle_region = ctx.declare_event_output<int>("region");
  handle_node = ctx.declare_event_output<int>("node"); //dnn region
  handle_num_b = ctx.declare_event_output<int>("num_b");
  handle_num_l = ctx.declare_event_output<int>("num_l");
  handle_num_j = ctx.declare_event_output<int>("num_j");
  handle_kfold = ctx.declare_event_output<int>("kfold"); 
}


bool InvReconstruction::process(Event& event){

  event.weight = event.get(handle_weight);

  //BTag scale factors
  sf_btagging->process(event);

  //Sanity check for lepton sf
  if (is_mc) {
    int leptons = (*event.electrons).size() + (*event.muons).size();
    if(leptons==0&&HasGenLepton(event)) h_missed->fill(event);
    if(leptons!=0&&(!HasGenLepton(event))) h_fake->fill(event);
  }

  //Global cuts
  h_base->fill(event);
  if (event.met->pt() < 170) return false;
  h_met->fill(event);
  if (DeltaPhi(event) < 0.5) return false;
  h_delta->fill(event);
  if (event.weight > 10) return false;
  h_weight->fill(event);

  //cutflow of signal region
  int lep = (*event.electrons).size() + (*event.muons).size();
  if (lep==0) h_veto->fill(event);
  if (lep==0 && !(s_bjet_none->passes(event))) h_btag->fill(event);
   
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

  int btags = 0;
  if (s_bjet_one->passes(event)) btags = 1;
  if (s_bjet_two->passes(event)) btags = 2;
  
  int jets = 5;
  if (s_njet_six->passes(event)) jets = 6;
  if (s_njet_sevenplus->passes(event)) jets = 7;

  event.set(handle_HT,HT);
  event.set(handle_met,event.met->pt());
  event.set(handle_num_l,lep);
  event.set(handle_num_b,btags);
  event.set(handle_num_j,jets);
  event.set(handle_delta_phi,DeltaPhi(event));
  event.set(handle_delta_eta,DeltaEta(event));
  event.set(handle_weight, event.weight);
  event.set(handle_kfold, get_random(event));
  event.set(handle_score, -1);
  event.set(handle_node, -1);

  return true;
}

bool InvReconstruction::AssignRegion(Event& event){
  int leptons = (*event.electrons).size() + (*event.muons).size();
  int btags = 0;
  if (s_bjet_one->passes(event)) btags = 1;
  if (s_bjet_two->passes(event)) btags = 2;
  
  int jets = 5;
  if (s_njet_sixplus->passes(event)) jets = 6;

  array<int,3> lbj = {leptons,btags,jets};
  Region region = Region::Invalid;

  if(lbj==array<int,3>{0,2,6}) region = Region::SR_6J;    //Signal region
  if(lbj==array<int,3>{0,2,5}) region = Region::SR_5J; 
  if(lbj==array<int,3>{0,1,6}) region = Region::SR_1B_6J; 
  if(lbj==array<int,3>{0,1,5}) region = Region::SR_1B_5J;
  if(lbj==array<int,3>{0,0,6}) region = Region::IR_0B_6J; //Invisible control region
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
  vector<Jet> jets = *event.jets;
  if(jets.size() < 1) return -1;

  double phi_m = event.met->phi();
  double min_diff = M_PI;

  for (Jet jet: jets){
    if(jet.pt()<30) continue;
    double phi_j = jet.phi();
    double delta = abs(phi_m - phi_j);
    if(delta > M_PI){delta = 2*M_PI - delta;}
    if(delta < min_diff) min_diff = delta;
  }

  return min_diff;
}

double InvReconstruction::DeltaEta(Event& event){
  vector<Jet> jets = *event.jets;
  if(jets.size() < 2) return -1 ;
  double eta1 = jets[0].eta();
  double eta2 = jets[1].eta();
  return abs(eta1-eta2);
}

double InvReconstruction::get_random(Event &event) {
  // Initialize random generator with phi-dependend random seed
  double phi = event.jets->at(0).v4().phi();
  srand((int)(1000 * phi));
  return (rand() % 5) + 1;
}

UHH2_REGISTER_ANALYSIS_MODULE(InvReconstruction)