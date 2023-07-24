#define CLASS HiggsReconstructor

#include <cmath>

#include "UHH2/AZH/include/HiggsReco.h"

using namespace std;
using namespace uhh2;

CLASS::HiggsReconstructor(Context & ctx) {
    reconstruction_hypothesis_builder.reset(new ReconstructionHypothesisBuilder(ctx));
    handle_ttbar_reco_hypotheses_fixedb = ctx.get_handle<vector<ReconstructionHypothesis>>("ttbar_reco_hypothese_fixedb");
}

bool CLASS::process(Event & event) {
  mt_A = 0;
  mt_H = 0;
  m_H = 0;

  bool has_five_jets = event.jets->size() == 5;

  LorentzVector H = {0,0,0,0};
  has_five_jets ? H = ProcessFiveJets(event) : H = ProcessMoreJets(event);

  const double met = event.met->pt();
  const double phi = event.met->phi();
  LorentzVector Z = LorentzVector(met,0,phi,sqrt(pow(z_mass,2)+pow(met,2)));
  LorentzVector Ht = Projection(H); 

  m_H = H.M();
  mt_A = (Z+Ht).M();

  return true;
}

LorentzVector CLASS::ProcessFiveJets(Event & event){
    LorentzVector H = {0,0,0,0};
    for (Jet jet : *event.jets) H += jet.v4();
    mt_H = -1; //alternative?
    return H;
}

LorentzVector CLASS::ProcessMoreJets(Event & event){
    reconstruction_hypothesis_builder->BuildRecoHypotheses(event);
    vector<ReconstructionHypothesis> hypotheses = event.get(handle_ttbar_reco_hypotheses_fixedb);
    ReconstructionHypothesis best_hypothesis = BestHypothesis(hypotheses);

    LorentzVector H_v4 = best_hypothesis.get_H_v4();
    LorentzVector t_v4 = best_hypothesis.get_t_v4();
    LorentzVector tbar_v4 = best_hypothesis.get_tbar_v4();
    
    mt_H = (Projection(t_v4)+Projection(tbar_v4)).M();
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

LorentzVector CLASS::Projection(LorentzVector v){
  const double pt = v.pt();
  const double phi = v.phi();
  const double m = v.M();  
  return LorentzVector(pt,0,phi,sqrt(pow(m,2)+pow(pt,2)));
}

void CLASS::TransMassAtGen(Event & event){
  LorentzVector h = LorentzVector();
  LorentzVector z = LorentzVector();
    for(GenParticle P : *event.genparticles){
      if (P.pdgId() == 23) z = LorentzVector(P.pt(),P.eta(),P.phi(),P.energy());
      else if(P.pdgId() == 35) h = LorentzVector(P.pt(),P.eta(),P.phi(),P.energy());
    }

    LorentzVector h_trans = Projection(h);
    LorentzVector z_trans = Projection(z);

    //gen_mass_H = h.M();
    //gen_mt_A = (h_trans+z_trans).M();
    //gen_met = z.pt();
}

double CLASS::GetTransMassA(){
  return mt_A;
}

double CLASS::GetTransMassH(){
  return mt_H;
}

double CLASS::GetMassH(){
  return m_H;
}

double CLASS::GetTransMassDiff(){
  return mt_A - mt_H;
}

double CLASS::GetMassDiff(){
  return mt_A - m_H;
}
