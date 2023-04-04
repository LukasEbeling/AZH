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


AtoZHRecoHists::AtoZHRecoHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // Handles
  handle_z_mass_reco = ctx.declare_event_output<double>("z_mass_reco");
  handle_z_pt_reco = ctx.declare_event_output<double>("z_pt_reco");
  handle_z_mass_gen_matched = ctx.declare_event_output<double>("z_mass_gen_matched");
  handle_z_pt_gen_matched = ctx.declare_event_output<double>("z_pt_gen_matched");
  handle_h_mass_reco_chi_sq = ctx.declare_event_output<double>("h_mass_reco_chi_sq");
  handle_h_mass_reco_six_jets = ctx.declare_event_output<double>("h_mass_reco_six_jets");
  handle_h_mass_gen_six_jets = ctx.declare_event_output<double>("h_mass_gen_six_jets");
  handle_a_mass = ctx.declare_event_output<double>("a_mass");
  handle_a_minus_h_mass = ctx.declare_event_output<double>("a_minus_h_mass");
  handle_tbar_mass_reco_chi_sq = ctx.declare_event_output<double>("tbar_mass_reco_chi_sq");

  // Histograms
  book<TH1F>("Z_Mass_Reco", "Z_{Mass}", 15, 80, 100);  
  book<TH1F>("Z_Pt_Reco", "Z_{Pt}", 30, 0, 700);  
  book<TH1F>("Z_Mass_Gen_Matched", "Z_{Mass, GenMatched}", 40, 80, 100);  
  book<TH1F>("Z_Pt_Gen_Matched", "Z_{Pt, GenMatched}", 30, 0, 700);  
  book<TH1F>("H_Mass_Reco_chi_sq", "H_{Mass, Chi2}", 30, 100, 1000); 
  book<TH1F>("H_Mass_Reco_six_jets", "H_{Mass, SixJets}", 30, 250, 900); 
  book<TH1F>("H_Mass_Gen_six_jets", "H_{Mass, GenSixJets}", 30, 250, 900); 
  book<TH1F>("TBar_Mass_Reco_chi_sq", "H_{Mass, Chi2}", 30, 0, 400); 
  book<TH1F>("A_Mass", "A_{Mass}", 30, 350, 1100); 
  book<TH1F>("deltaM", "deltaM", 22, 0, 700); 
}


void AtoZHRecoHists::fill(const Event &event){
  double weight = event.weight;

  double z_mass_reco = event.get(handle_z_mass_reco);
  double z_mass_gen_matched = event.get(handle_z_mass_gen_matched);
  double z_pt_reco = event.get(handle_z_pt_reco);
  double z_pt_gen_matched = event.get(handle_z_pt_gen_matched);

  double h_mass_reco_chi_sq = event.get(handle_h_mass_reco_chi_sq);
  double h_mass_reco_six_jets = event.get(handle_h_mass_reco_six_jets);
  double h_mass_gen_six_jets = event.get(handle_h_mass_gen_six_jets);

  double a_mass = event.get(handle_a_mass);
  double a_minus_h_mass = event.get(handle_a_minus_h_mass);

  double tbar_mass_reco_chi_sq = event.get(handle_tbar_mass_reco_chi_sq);

  hist("Z_Mass_Reco")->Fill(z_mass_reco, weight);
  hist("Z_Pt_Reco")->Fill(z_pt_reco, weight);
  hist("Z_Mass_Gen_Matched")->Fill(z_mass_gen_matched, weight);
  hist("Z_Pt_Gen_Matched")->Fill(z_pt_gen_matched, weight);
  hist("H_Mass_Reco_chi_sq")->Fill(h_mass_reco_chi_sq, weight);
  hist("H_Mass_Reco_six_jets")->Fill(h_mass_reco_six_jets, weight);
  hist("H_Mass_Gen_six_jets")->Fill(h_mass_gen_six_jets, weight);
  hist("A_Mass")->Fill(a_mass, weight);
  hist("deltaM")->Fill(a_minus_h_mass, weight);
  hist("TBar_Mass_Reco_chi_sq")->Fill(tbar_mass_reco_chi_sq, weight);
}


AtoZHRecoHists::~AtoZHRecoHists(){}


AtoZHZMassWindowScanHists::AtoZHZMassWindowScanHists(Context & ctx, const string & dirname): Hists(ctx, dirname) {
  // Handles
  handle_z_mass_reco = ctx.declare_event_output<double>("z_mass_reco");
  // Histograms
  book<TH1F>("Z_Mass_Window", "Z_{Mass,Window}", 12, 1, 13);  
}


void AtoZHZMassWindowScanHists::fill(const Event & event){
  double weight = event.weight;
  double z_mass = event.get(handle_z_mass_reco);
  double delta = abs(z_mass - 91.188);

  if (delta < 1) { hist("Z_Mass_Window")->Fill(1, weight); }
  if (delta < 2) { hist("Z_Mass_Window")->Fill(2, weight); }
  if (delta < 3) { hist("Z_Mass_Window")->Fill(3, weight); }
  if (delta < 4) { hist("Z_Mass_Window")->Fill(4, weight); }
  if (delta < 5) { hist("Z_Mass_Window")->Fill(5, weight); }
  if (delta < 6) { hist("Z_Mass_Window")->Fill(6, weight); }
  if (delta < 7) { hist("Z_Mass_Window")->Fill(7, weight); }
  if (delta < 8) { hist("Z_Mass_Window")->Fill(8, weight); }
  if (delta < 9) { hist("Z_Mass_Window")->Fill(9, weight); }
  if (delta < 10) { hist("Z_Mass_Window")->Fill(10, weight); }
  if (delta < 11) { hist("Z_Mass_Window")->Fill(11, weight); }
  if (delta < 12) { hist("Z_Mass_Window")->Fill(12, weight); }
}


AtoZHZMassWindowScanHists::~AtoZHZMassWindowScanHists(){}


AtoZHHists::AtoZHHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // selection yields
  book<TH1F>("N_Events", "N_{Events}", 1, -1, 1);  

  // jets
  book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);  
  book<TH1F>("N_PU", "N_{PU}", 100, 0, 100);  
  book<TH1F>("eta_jet1", "#eta^{jet 1}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet2", "#eta^{jet 2}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet3", "#eta^{jet 3}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet4", "#eta^{jet 4}", 40, -2.5, 2.5);
  book<TH1F>("phi_jet1", "#phi^{jet 1}", 40, -2.5, 2.5);
  book<TH1F>("pt_jet1", "p_{T}^{jet 1}", 100, 10, 500);
  book<TH1F>("pt_jet2", "p_{T}^{jet 2}", 100, 10, 500);
  book<TH1F>("pt_jet3", "p_{T}^{jet 3}", 100, 10, 500);
  book<TH1F>("pt_jet4", "p_{T}^{jet 4}", 100, 10, 500);

  book<TH1F>("EMcharged_jet1", "EMcharged_jet1", 100,0.0,1.0);
  book<TH1F>("EMneutral_jet1", "EMneutral_jet1", 100,0.0,1.0);
  book<TH1F>("HADcharged_jet1", "HADcharged_jet1", 100,0.0,1.0);
  book<TH1F>("HADneutral_jet1", "HADneutral_jet1", 100,0.0,1.0);

  book<TH2D>("EMcharged_vs_eta_jet1","EMcharged vs #eta; #eta; EMcharged",100,-6,6,100,0.0,1.0);   
  book<TH2D>("EMneutral_vs_eta_jet1","EMneutral vs #eta; #eta; EMneutral",100,-6,6,100,0.0,1.0);   
  book<TH2D>("HADcharged_vs_eta_jet1","HADcharged vs #eta; #eta; HADcharged",100,-6,6,100,0.0,1.0);   
  book<TH2D>("HADneutral_vs_eta_jet1","HADneutral vs #eta; #eta; HADneutral",100,-6,6,100,0.0,1.0);   
  book<TH2D>("EMcharged_vs_PU_jet1","EMcharged vs PU; PU; EMcharged",100,0,100,100,0.0,1.0);   
  book<TH2D>("EMneutral_vs_PU_jet1","EMneutral vs PU; PU; EMneutral",100,0,100,100,0.0,1.0);   
  book<TH2D>("HADcharged_vs_PU_jet1","HADcharged vs PU; PU; HADcharged",100,0,100,100,0.0,1.0);   
  book<TH2D>("HADneutral_vs_PU_jet1","HADneutral vs PU; PU; HADneutral",100,0,100,100,0.0,1.0);   

  book<TH1F>("n_bjets_loose", "N_{b jets, loose}", 10, 0, 10);
  book<TH1F>("n_bjets_medium", "N_{b jets, medium}", 10, 0, 10);
  book<TH1F>("n_bjets_tight", "N_{b jets, tight}", 10, 0, 10);

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 40, 0, 200);
  book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  book<TH1F>("N_electrons", "N^{#electrons}", 10, 0, 10);
  book<TH1F>("pt_electrons", "p_{T}^{#e} [GeV/c]", 40, 0, 200);
  book<TH1F>("eta_electrons", "#eta^{#e}", 40, -2.1, 2.1);
  book<TH1F>("reliso_electrons", "#e rel. Iso", 40, 0, 0.5);
  book<TH1F>("missing_pt","Missing pt", 200, 0, 1000);

  // primary vertices
  book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);
}


void AtoZHHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  hist("N_Events")->Fill(0., weight);

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

  hist("missing_pt")->Fill(event.met->pt(), weight);
  hist("n_bjets_loose")->Fill(n_btag_jets_loose, weight);
  hist("n_bjets_medium")->Fill(n_btag_jets_medium, weight);
  hist("n_bjets_tight")->Fill(n_btag_jets_tight, weight);

  int n_jets = jets->size();
  hist("N_jets")->Fill(n_jets, weight);
  if(!event.isRealData) hist("N_PU")->Fill(event.genInfo->pileup_TrueNumInteractions(), weight);

  if(n_jets>=1){
    hist("eta_jet1")->Fill(jets->at(0).eta(), weight);
    hist("phi_jet1")->Fill(jets->at(0).phi(), weight);
    hist("pt_jet1")->Fill(jets->at(0).pt(), weight);
    hist("EMcharged_jet1")->Fill(jets->at(0).chargedEmEnergyFraction(), weight);
    hist("EMneutral_jet1")->Fill(jets->at(0).neutralEmEnergyFraction(), weight);
    hist("HADcharged_jet1")->Fill(jets->at(0).chargedHadronEnergyFraction(), weight);
    hist("HADneutral_jet1")->Fill(jets->at(0).neutralHadronEnergyFraction(), weight);
    
    ((TH2D*)hist("EMcharged_vs_eta_jet1"))->Fill(jets->at(0).eta(), jets->at(0).chargedEmEnergyFraction(), weight);
    ((TH2D*)hist("EMneutral_vs_eta_jet1"))->Fill(jets->at(0).eta(), jets->at(0).neutralEmEnergyFraction(), weight);
    ((TH2D*)hist("HADcharged_vs_eta_jet1"))->Fill(jets->at(0).eta(), jets->at(0).chargedHadronEnergyFraction(), weight);
    ((TH2D*)hist("HADneutral_vs_eta_jet1"))->Fill(jets->at(0).eta(), jets->at(0).neutralHadronEnergyFraction(), weight);
    if(!event.isRealData){
      ((TH2D*)hist("EMcharged_vs_PU_jet1"))->Fill(event.genInfo->pileup_TrueNumInteractions(), jets->at(0).chargedEmEnergyFraction(), weight);
      ((TH2D*)hist("EMneutral_vs_PU_jet1"))->Fill(event.genInfo->pileup_TrueNumInteractions(), jets->at(0).neutralEmEnergyFraction(), weight);
      ((TH2D*)hist("HADcharged_vs_PU_jet1"))->Fill(event.genInfo->pileup_TrueNumInteractions(), jets->at(0).chargedHadronEnergyFraction(), weight);
      ((TH2D*)hist("HADneutral_vs_PU_jet1"))->Fill(event.genInfo->pileup_TrueNumInteractions(), jets->at(0).neutralHadronEnergyFraction(), weight);
    }
  }
  if(n_jets>=2){
    hist("eta_jet2")->Fill(jets->at(1).eta(), weight);
    hist("pt_jet2")->Fill(jets->at(1).pt(), weight);
  }
  if(n_jets>=3){
    hist("eta_jet3")->Fill(jets->at(2).eta(), weight);
    hist("pt_jet3")->Fill(jets->at(2).pt(), weight);
  }
  if(n_jets>=4){
    hist("eta_jet4")->Fill(jets->at(3).eta(), weight);
    hist("pt_jet4")->Fill(jets->at(3).pt(), weight);
  }

  int n_muons = event.muons->size();
  hist("N_mu")->Fill(n_muons, weight);
  for (const Muon & thismu : *event.muons){
    hist("pt_mu")->Fill(thismu.pt(), weight);
    hist("eta_mu")->Fill(thismu.eta(), weight);
    hist("reliso_mu")->Fill(thismu.relIso(), weight);
  }

  int n_electrons = event.electrons->size();
  hist("N_electrons")->Fill(n_electrons, weight);
  for (const Electron & ele : *event.electrons){
    hist("pt_electrons")->Fill(ele.pt(), weight);
    hist("eta_electrons")->Fill(ele.eta(), weight);
    hist("reliso_electrons")->Fill(ele.relIso(), weight);
  }

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);
}

AtoZHHists::~AtoZHHists(){}


StandardModelBRHists::StandardModelBRHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // Histograms
  book<TH1F>("ee_mumu", "ee / mumu", 2, 0.5, 2.5);  
  book<TH1F>("ee_mumu_weighted", "ee / mumu", 2, 0.5, 2.5);  
  book<TH1F>("W_lep_had", "W lep / had", 2, 0.5, 2.5);  
  book<TH1F>("W_lep_had_weighted", "W lep / had", 2, 0.5, 2.5);  
  book<TH1F>("N_gen_Z", "N_{gen}(Z)", 1, 0.5, 2.5);  
  book<TH1F>("N_fully_leptonic", "N fully leptonic", 1, 0.5, 2.5);  
  book<TH1F>("N_gen_Wlep", "N_{gen}(W_{lep})", 1, 0.5, 2.5);  
}

void StandardModelBRHists::fill(const Event & event){
  double weight = event.weight;
  const vector<GenParticle> *genparticles = event.genparticles;

  // ee mumu
  for (unsigned i=0; i<genparticles->size(); i++){
    const GenParticle gp = genparticles->at(i);
    if ( !is_emu(gp) ) continue;
    const GenParticle m1 = *genparticles->at(i).mother(genparticles, 1);
    if ( !is_Z(m1) ) continue;
    hist("N_gen_Z")->Fill(1, weight);
    if ( is_Muon(gp) ) {
      hist("ee_mumu")->Fill(1, 1.);
      hist("ee_mumu_weighted")->Fill(1, weight);
    }
    else {
      hist("ee_mumu")->Fill(2, 1.);
      hist("ee_mumu_weighted")->Fill(2, weight);
    }
  }

  // W
  for (unsigned i=0; i<genparticles->size(); i++){
    const GenParticle gp = genparticles->at(i);
    if ( !is_W(gp) ) continue;
    const GenParticle d1 = *genparticles->at(i).daughter(genparticles, 1);
    const GenParticle d2 = *genparticles->at(i).daughter(genparticles, 2);

    if ( is_charged_lepton(d1) || is_charged_lepton(d2)) {
      hist("N_gen_Wlep")->Fill(1, weight);
      hist("W_lep_had_weighted")->Fill(1, weight);
    }
    if ( is_emu(d1) && is_emu(d2)) {
      hist("N_fully_leptonic")->Fill(1, weight);
    }
    if ( is_light_quark(d1) || is_b_bbar(d2) ) {
      hist("W_lep_had_weighted")->Fill(2, weight);
    }
  }
}

StandardModelBRHists::~StandardModelBRHists(){}


TriggerHists::TriggerHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // Handles
  handle_channel = ctx.get_handle<int>("channel");
  handle_muons_tight = ctx.get_handle<vector<Muon>>("muons_tight");
  handle_electrons_tight = ctx.get_handle<vector<Electron>>("electrons_tight");

  // Histograms
  book<TH1F>("diEleChannel_1", "p_{T}(e_1)", 25, 0, 800);
  book<TH1F>("diEleChannel_2", "p_{T}(e_2)", 25, 0, 800);
  book<TH1F>("diMuChannel_1", "p_{T}(mu_1)", 25, 0, 800);
  book<TH1F>("diMuChannel_2", "p_{T}(mu_2)", 25, 0, 800);
  book<TH1F>("EleMuChannel_ele", "p_{T}(e)", 25, 0, 800);
  book<TH1F>("EleMuChannel_mu", "p_{T}(mu)", 25, 0, 800);
}

void TriggerHists::fill(const Event &event){
  double weight = event.weight;
  vector<Muon> muons_tight = event.get(handle_muons_tight);
  vector<Electron> electrons_tight = event.get(handle_electrons_tight);
  int x = event.get(handle_channel);
  Channel channel = static_cast<Channel>(x);

  if ( channel == Channel::diElectron ) {
    hist("diEleChannel_1")->Fill(electrons_tight.at(0).pt(), weight);
    hist("diEleChannel_2")->Fill(electrons_tight.at(1).pt(), weight);
  }
  if ( channel == Channel::diMuon ) {
    hist("diMuChannel_1")->Fill(muons_tight.at(0).pt(), weight);
    hist("diMuChannel_2")->Fill(muons_tight.at(1).pt(), weight);
  }
  if ( channel == Channel::ElectronMuon ) {
    hist("EleMuChannel_ele")->Fill(electrons_tight.at(0).pt(), weight);
    hist("EleMuChannel_mu")->Fill(muons_tight.at(0).pt(), weight);
  }
}

TriggerHists::~TriggerHists(){}

