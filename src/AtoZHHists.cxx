#include "UHH2/AZH/include/AtoZHHists.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/AZH/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <math.h>

using namespace std;
using namespace uhh2;
using namespace pdgIdUtils;


HistSet::HistSet(Context & ctx, const string & dirname): Hists(ctx, dirname){

  book<TH1F>("N_Events", "N_{Events}", 1, -1, 1);  
  book<TH1F>("missing_pt","Missing pt", int(1500/20), 0, 1500);
  book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);
  book<TH1F>("N_lep", "N_{lep}", 5, 0, 5);
  book<TH1F>("n_bjets_tight", "N_{b jets, tight}", 5, 0, 5);
  book<TH1F>("n_bjets_loose", "N_{b jets, loose}", 5, 0, 5);
  book<TH1F>("n_bjets_medium", "N_{b jets, medium}", 5, 0, 5);
  book<TH1F>("delta_phi", "delta_phi",10,0,1*M_PI);    
}


void HistSet::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  hist("N_Events")->Fill(0., weight);
  hist("missing_pt")->Fill(event.met->pt(), weight);


  std::vector<Jet>* jets = event.jets;

  const BTag::algo btag_algo = BTag::DEEPJET;
  const JetId btag_id_loose = BTag(btag_algo, BTag::WP_LOOSE);
  const JetId btag_id_medium = BTag(btag_algo, BTag::WP_MEDIUM);
  const JetId btag_id_tight = BTag(btag_algo, BTag::WP_TIGHT);

  int n_btag_jets_loose = 0;
  int n_btag_jets_medium = 0;
  int n_btag_jets_tight = 0;
  for(const Jet & jet : *event.jets) {
    if(btag_id_loose(jet, event)) ++n_btag_jets_loose;
    if(btag_id_medium(jet, event)) ++n_btag_jets_medium;
    if(btag_id_tight(jet, event)) ++n_btag_jets_tight;
  }

  hist("N_jets")->Fill(jets->size(), weight);
  hist("n_bjets_loose")->Fill(n_btag_jets_loose, weight);
  hist("n_bjets_medium")->Fill(n_btag_jets_medium, weight);
  hist("n_bjets_tight")->Fill(n_btag_jets_tight, weight);
  hist("N_lep")->Fill((*event.electrons).size() + (*event.muons).size(),weight);
  hist("delta_phi")->Fill(DeltaPhi(event),weight);
}

HistSet::~HistSet(){}

// calculate Delta Phi (min angle between met and jets)
double HistSet::DeltaPhi(const Event& event){
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

// Set of Histogram only with total number of events
SimpleHist::SimpleHist(Context & ctx, const string & dirname): Hists(ctx, dirname){

  book<TH1F>("N_Events", "N_{Events}", 1, -1, 1);
}

void SimpleHist::fill(const Event & event){
  double weight = event.weight;
  hist("N_Events")->Fill(0., weight);
}

SimpleHist::~SimpleHist(){}
