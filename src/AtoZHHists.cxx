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
  //hist("n_bjets_loose")->Fill(n_btag_jets_loose, weight);
  //hist("n_bjets_medium")->Fill(n_btag_jets_medium, weight);
  //hist("n_bjets_tight")->Fill(n_btag_jets_tight, weight);
}

AtoZHHists::~AtoZHHists(){}