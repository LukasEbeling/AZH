#include <iostream>
#include <cstdlib>
#include <cmath>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/PSWeights.h"

#include "UHH2/AZH/include/JetCleaner.h"
#include "UHH2/AZH/include/NormalisationTools.h"
#include "UHH2/AZH/include/RochesterCorrections.h"
#include "UHH2/AZH/include/ScaleFactors.h"
#include "UHH2/AZH/include/Utils.h"
#include "UHH2/AZH/include/HiggsReco.h"


using namespace std;
using namespace uhh2;


// Prefix h_ denotes histogram
// Prefix s_ denotes selection
// Prefix handle_ denotes handle


class InvTopology: public AnalysisModule {
  public:
    explicit InvTopology(Context & ctx);
    virtual bool process(Event & event) override;
    double DeltaPhi(Event & event);

  private:
    bool is_mc;
    Year year;
    const JetId jet_id = AndId<Jet>(PtEtaCut(30, 2.4), JetPFID(JetPFID::WP_TIGHT_LEPVETO_CHS));
    MuonId muonId_loose;
    ElectronId electronId_loose;  

    // Handles
    uhh2::Event::Handle<double> handle_event_weight;
    uhh2::Event::Handle<double> handle_origin_weight;
    uhh2::Event::Handle<double> handle_delta_phi;
    uhh2::Event::Handle<double> handle_met;
    uhh2::Event::Handle<double> handle_A_mt;
    uhh2::Event::Handle<double> handle_H_m;
    uhh2::Event::Handle<int> handle_num_jets;
    uhh2::Event::Handle<int> handle_num_lep;
    uhh2::Event::Handle<int> handle_num_btag;


    // Rochester corrections and leptons ID+ISO collections
    std::unique_ptr<AnalysisModule> rochester;

    // Selection
    std::unique_ptr<CommonModules> common_modules;
    std::unique_ptr<AnalysisModule> clnr_jetpuid;

    // Studies on jet number
    std::unique_ptr<Selection> s_njet_none;
    std::unique_ptr<Selection> s_njet_one;
    std::unique_ptr<Selection> s_njet_two;
    std::unique_ptr<Selection> s_njet_three;
    std::unique_ptr<Selection> s_njet_four;
    std::unique_ptr<Selection> s_njet_five;
    std::unique_ptr<Selection> s_njet_six;
    std::unique_ptr<Selection> s_njet_seven;
    std::unique_ptr<Selection> s_bjet_none;
    std::unique_ptr<Selection> s_bjet_one;
    std::unique_ptr<Selection> s_bjet_two;

    // Event Weighting
    unique_ptr<HEMSelection> sel_hem;
    unique_ptr<MCLumiWeight> sf_lumi;
    unique_ptr<MCPileupReweight> sf_pileup;
    unique_ptr<AnalysisModule> sf_l1prefiring;
    unique_ptr<AnalysisModule> sf_vjets;
    unique_ptr<AnalysisModule> sf_mtop;
    unique_ptr<AnalysisModule> sf_puid;
    unique_ptr<MCScaleVariation> sf_QCDScaleVar;
    unique_ptr<PSWeights> ps_weights;
    unique_ptr<AnalysisModule> pdf_weights;

    std::unique_ptr<HiggsReconstructor> higgs_reconstructor;
}; 


InvTopology::InvTopology(Context & ctx){
  year = extract_year(ctx);
  is_mc = ctx.get("dataset_type") == "MC";
  higgs_reconstructor.reset(new HiggsReconstructor(ctx));

  // Cleaners
  const JetId bmedium = BTag(BTag::DEEPJET, BTag::WP_MEDIUM);
  const JetId btight = BTag(BTag::DEEPJET, BTag::WP_TIGHT);

  // MuonId
  Muon::Selector muonIdTag = Muon::CutBasedIdLoose;
  Muon::Selector muonIsoTag = Muon::PFIsoLoose;
  muonId_loose = AndId<Muon>(PtEtaCut(20, 2.4), MuonID(muonIdTag), MuonID(muonIsoTag));

  // ElectronId
  Electron::tag eleTag = Electron::mvaEleID_Fall17_iso_V2_wp90;
  electronId_loose = AndId<Electron>(ElectronEtaWindowId(), PtEtaSCCut(20, 2.4), ElectronTagID(eleTag));

  // Selections
  s_njet_none.reset(new NJetSelection(0,0));
  s_njet_one.reset(new NJetSelection(1,1));
  s_njet_two.reset(new NJetSelection(2,2));
  s_njet_three.reset(new NJetSelection(3,3));
  s_njet_four.reset(new NJetSelection(4,4));
  s_njet_five.reset(new NJetSelection(5,5));
  s_njet_six.reset(new NJetSelection(6,6));
  s_njet_seven.reset(new NJetSelection(7,7));
  s_bjet_none.reset(new NJetSelection(0,0,bmedium));
  s_bjet_one.reset(new NJetSelection(1,1,bmedium));
  s_bjet_two.reset(new NJetSelection(2,2,bmedium));

  // Event Weighting
  sel_hem.reset(new HEMSelection(ctx));
  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_pileup.reset(new MCPileupReweight(ctx));
  sf_l1prefiring.reset(new L1PrefiringWeight(ctx));
  sf_vjets.reset(new VJetsReweighting(ctx));
  sf_mtop.reset(new TopPtReweighting(ctx, string2bool(ctx.get("apply_TopPtReweighting"))));
  sf_QCDScaleVar.reset(new MCScaleVariation(ctx));
  sf_puid.reset(new PUIDScaleFactors(ctx));
  pdf_weights.reset(new PDFWeightHandleProducer(ctx)); 
  ps_weights.reset(new PSWeights(ctx));

  // Handles
  handle_event_weight = ctx.declare_event_output<double>("event_weight");
  handle_origin_weight = ctx.declare_event_output<double>("origin_weight");
  handle_delta_phi = ctx.declare_event_output<double>("delta_phi");
  handle_met = ctx.declare_event_output<double>("met");
  handle_A_mt = ctx.declare_event_output<double>("mt_A");
  handle_H_m = ctx.declare_event_output<double>("m_H");
  handle_num_jets = ctx.declare_event_output<int>("num_jets");
  handle_num_lep = ctx.declare_event_output<int>("num_leps");
  handle_num_btag = ctx.declare_event_output<int>("num_btags");


  // Common
  common_modules.reset(new CommonModules());
  common_modules-> disable_mclumiweight(); //done manually 
  common_modules-> disable_mcpileupreweight(); //done manually 
  common_modules->set_jet_id(jet_id);
  common_modules->switch_jetlepcleaner(true);
  common_modules->switch_jetPtSorter();
  //common_modules-> disable_lumisel();
  //common_modules-> disable_metfilters();
  //common_modules-> disable_pvfilter();
  //common_modules->disable_jec();
  //common_modules->disable_jersmear();
  //common_modules->disable_jetpfidfilter();
  common_modules->init(ctx);

  rochester.reset(new RochesterCorrections(ctx));
  clnr_jetpuid.reset(new JetCleaner(ctx, JetPUID(JetPUID::WP_TIGHT)));
}

bool InvTopology::process(Event & event) {
  
  event.set(handle_origin_weight,event.weight);

  // Apply Jet corrections
  bool passes_common = common_modules->process(event);
  if (!passes_common) { return false; }

  // Apply rochester corrections before ID'ing muons
  rochester->process(event);

  //Lepton Cleaning
  clean_collection(*event.muons, event, muonId_loose);
  clean_collection(*event.electrons, event, electronId_loose);  

  // Jet PU ID, to be applied on JEC corrected jets
  // JECs in common_modules, hence here PU ID
  bool pileup_cleaned = clnr_jetpuid->process(event);
  if (!pileup_cleaned) {return false;}

  //HEM issue
  if(sel_hem->passes(event)) {
    if(event.isRealData) return false;
    else event.weight *= (1. - sel_hem->GetAffectedLumiFraction());
  }  

  // Event Weighting
  sf_lumi->process(event);
  sf_pileup->process(event);
  sf_l1prefiring->process(event);
  sf_vjets->process(event);
  sf_mtop->process(event); 
  sf_QCDScaleVar->process(event);
  sf_puid->process(event);
  pdf_weights->process(event);
  ps_weights->process(event);

  //number of jets
  int num_jets = 8; //eight or more
  if (s_njet_none->passes(event)) {num_jets = 0;}
  if (s_njet_one->passes(event)) {num_jets = 1;}
  if (s_njet_two->passes(event)) {num_jets = 2;}
  if (s_njet_three->passes(event)) {num_jets = 3;}
  if (s_njet_four->passes(event)) {num_jets = 4;}
  if (s_njet_five->passes(event)) {num_jets = 5;}
  if (s_njet_six->passes(event)) {num_jets = 6;}
  if (s_njet_seven->passes(event)) {num_jets = 7;}
  event.set(handle_num_jets, num_jets);

  //number of leptons
  int num_lep = (*event.electrons).size() + (*event.muons).size();
  event.set(handle_num_lep, num_lep);

  //number of b-tags
  int btags = 3; //three or more
  if (s_bjet_none->passes(event)) btags = 0;
  if (s_bjet_one->passes(event)) btags = 1;
  if (s_bjet_two->passes(event)) btags = 2;
  event.set(handle_num_btag,btags);

  //missing ET
  event.set(handle_met,event.met->pt());

  //delta phi min
  event.set(handle_delta_phi,DeltaPhi(event));

  //Reconstruct A and H
  double mt_A = -1;
  double m_H = -1;

  if (num_jets != 0) {
    higgs_reconstructor->process(event);
    mt_A = higgs_reconstructor->GetTransMassA();
    m_H = higgs_reconstructor->GetMassH();
  }
  
  event.set(handle_A_mt,mt_A);
  event.set(handle_H_m,m_H);

  event.set(handle_event_weight, event.weight);

  return true;
}

double InvTopology::DeltaPhi(Event& event){
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

UHH2_REGISTER_ANALYSIS_MODULE(InvTopology)

