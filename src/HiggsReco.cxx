#define CLASS HiggsReconstructor

#include <cmath>

#include "UHH2/AZH/include/HiggsReco.h"

using namespace std;
using namespace uhh2;

CLASS::HiggsReconstructor(Context & ctx) {
    reconstruction_hypothesis_builder.reset(new ReconstructionHypothesisBuilder(ctx));
    handle_ttbar_reco_hypotheses_fixedb = ctx.get_handle<vector<ReconstructionHypothesis>>("ttbar_reco_hypothese_fixedb");
    handle_H_m = ctx.declare_event_output<double>("H_m");
    handle_met = ctx.declare_event_output<double>("MET");
    handle_A_mt = ctx.declare_event_output<double>("A_mt");
    handle_H_mt = ctx.declare_event_output<double>("H_mt");
    handle_W_m = ctx.declare_event_output<double>("W_m");
    handle_gen_A_mt = ctx.declare_event_output<double>("A_mt_gen");
}

bool CLASS::process(Event & event) {

    event.set(handle_H_m,0);
    event.set(handle_met,0);
    event.set(handle_A_mt,0);
    //get number of jets, replace later by regions
    bool has_five_jets = event.jets->size() == 5;

    LorentzVector H = {0,0,0,0};
    has_five_jets ? H = ProcessFiveJets(event) : H = ProcessMoreJets(event);

    const double met = event.met->pt();
    const double phi = event.met->phi();
    LorentzVector Z = LorentzVector(met,0,phi,sqrt(pow(z_mass,2)+pow(met,2)));
    LorentzVector Ht = Projection(H); 

    event.set(handle_A_mt,(Z+Ht).M());    
    event.set(handle_met,met);
    event.set(handle_H_m,H.M());
    event.set(handle_gen_A_mt,TransMassAtGen(event));

    return true;
}

LorentzVector CLASS::ProcessFiveJets(Event & event){
    LorentzVector H = {0,0,0,0};
    for (Jet jet : *event.jets) H += Projection(jet.v4());
    event.set(handle_H_mt,-1); //alternative?
    event.set(handle_W_m,-1); //alternative?
    return H;
}

LorentzVector CLASS::ProcessMoreJets(Event & event){
    reconstruction_hypothesis_builder->BuildRecoHypotheses(event);
    vector<ReconstructionHypothesis> hypotheses = event.get(handle_ttbar_reco_hypotheses_fixedb);
    ReconstructionHypothesis best_hypothesis = BestHypothesis(hypotheses);
    
    event.set(handle_W_m,BestWMass(hypotheses));

    LorentzVector H_v4 = best_hypothesis.get_H_v4();
    LorentzVector t_v4 = best_hypothesis.get_t_v4();
    LorentzVector tbar_v4 = best_hypothesis.get_tbar_v4();
    const double mt = (Projection(t_v4)+Projection(tbar_v4)).M();
    event.set(handle_H_mt,mt);
    return H_v4;
}

ReconstructionHypothesis CLASS::BestHypothesis(vector<ReconstructionHypothesis> hypotheses) {
  ReconstructionHypothesis best_hypothesis = hypotheses.at(0);
  double lowest_chi_sq = ChiSquared(best_hypothesis);
  double chi_sq;
  for (unsigned i=1; i<hypotheses.size(); i++) {
    chi_sq = ChiSquared(hypotheses.at(i));
    if (chi_sq < lowest_chi_sq) {
      lowest_chi_sq = chi_sq;
      best_hypothesis = hypotheses.at(i);
    }
  }
  return best_hypothesis;
}

double CLASS::ChiSquared(ReconstructionHypothesis hypothesis){
  double m_t1 = hypothesis.get_t_m();
  double m_w1 = hypothesis.get_t_w_m();
  double m_t2 = hypothesis.get_tbar_m();
  double m_w2 = hypothesis.get_tbar_w_m();
  double chi_sq = pow((m_t1 - top_mass) / sigma_top, 2)
    + pow((m_t2 - top_mass) / sigma_top, 2)
    + pow((m_w1 - w_mass) / sigma_w, 2)
    + pow((m_w2 - w_mass) / sigma_w, 2);
  return chi_sq;
}

double CLASS::BestWMass(vector<ReconstructionHypothesis> hypotheses){
double best_w_mass = WMass(hypotheses.at(0));
for (ReconstructionHypothesis hyp : hypotheses){
  double w_mass = WMass(hyp);
  if (abs(w_mass-80.38) < abs(best_w_mass-80.38)) best_w_mass = w_mass;
}
return best_w_mass;
}

double CLASS::WMass(ReconstructionHypothesis hypothesis){
  LorentzVector q1 = hypothesis.get_t_q1_v4();
  LorentzVector q2 = hypothesis.get_t_q2_v4();
  LorentzVector q3 = hypothesis.get_tbar_q1_v4();
  LorentzVector q4 = hypothesis.get_tbar_q2_v4();
  double w1 = (q1+q2).M();
  double w2 = (q3+q4).M();
  return abs(w1-80.38) > abs(w2-80.38) ? w1 : w2;
}

LorentzVector CLASS::Projection(LorentzVector v){
  const double pt = v.pt();
  const double phi = v.phi();
  const double m = v.M();  
  return LorentzVector(pt,0,phi,sqrt(pow(m,2)+pow(pt,2)));
}


double CLASS::TransMassAtGen(Event & event){
  LorentzVector h = LorentzVector();
  LorentzVector z = LorentzVector();
    for(GenParticle P : *event.genparticles){
      if (P.pdgId() == 23) z = LorentzVector(P.pt(),P.eta(),P.phi(),P.energy());
      else if(P.pdgId() == 35) h = LorentzVector(P.pt(),P.eta(),P.phi(),P.energy());
    }
    h = Projection(h);
    z = Projection(z);
    return (h+z).M();
}

/*
double CLASS::TransMassAtGen(Event & event){
  LorentzVector H = {0,0,0,0};
  LorentzVector Z = {0,0,0,0};
  for(GenParticle P : *event.genparticles){
      if (P.pdgId() == 23) Z = LorentzVector(P.pt(),P.eta(),P.phi(),P.energy());
  }
  for(GenJet jet : *event.genjets) H+=jet.v4();
  //H = Projection(H);
  //Z = Projection(Z);
  return (H).M();
}
*/