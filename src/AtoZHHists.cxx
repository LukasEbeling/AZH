#include "UHH2/AZH/include/AtoZHHists.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/AZH/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace pdgIdUtils;


AtoZHHists::AtoZHHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  book<TH1F>("N_Events", "N_{Events}", 1, -1, 1);  
  book<TH1F>("missing_pt","Missing pt", 200, 0, 1000);
  book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);
  book<TH1F>("n_bjets_tight", "N_{b jets, tight}", 10, 0, 10);  
}


void AtoZHHists::fill(const Event & event){
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
}

AtoZHHists::~AtoZHHists(){}