#include <iostream>
#include <cstdlib>
#include <cmath>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/AZH/include/AtoZHHists.h"
#include "UHH2/AZH/include/JetCleaner.h"
#include "UHH2/AZH/include/MetFilters.h"
#include "UHH2/AZH/include/NormalisationTools.h"
#include "UHH2/AZH/include/RochesterCorrections.h"
#include "UHH2/AZH/include/ScaleFactors.h"
#include "UHH2/AZH/include/Triggers.h"
#include "UHH2/AZH/include/Utils.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"


using namespace std;
using namespace uhh2;


// Prefix h_ denotes histogram
// Prefix s_ denotes selection
// Prefix handle_ denotes handle


class InvPreselection: public AnalysisModule {
  public:
    explicit InvPreselection(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    bool is_mc;
    Year year;
    const JetId jet_id = AndId<Jet>(PtEtaCut(30, 2.4), JetPFID(JetPFID::WP_TIGHT_LEPVETO_CHS));
    MuonId muonId_loose;
    ElectronId electronId_loose;

    // Selection niherits from AnalysisModule as it uses passes() and process()
    std::unique_ptr<AnalysisModule> s_metFilters;

    // Handles
    uhh2::Event::Handle<double> handle_event_weight;
    uhh2::Event::Handle<int> handle_tight_b;

    // Histograms
    std::unique_ptr<Hists> h_unc_norm;
    std::unique_ptr<Hists> h_baseline;
    std::unique_ptr<Hists> h_no_leptons;
    std::unique_ptr<Hists> h_six_jets;
    //std::unique_ptr<Hists> h_bjet_none;
    //std::unique_ptr<Hists> h_bjet_one;
    std::unique_ptr<Hists> h_bjet_two;
    std::unique_ptr<Hists> h_met_cut;
    std::unique_ptr<Hists> h_delta_cut;


    // Histograms for BTagging efficiency measurements
    std::unique_ptr<BTagMCEfficiencyHists> h_btag_eff;

    // Rochester corrections and leptons ID+ISO collections
    std::unique_ptr<AnalysisModule> rochester;

    // Other Modules
    std::unique_ptr<CommonModules> common_modules;
    std::unique_ptr<AnalysisModule> clnr_jetpuid;
    std::unique_ptr<Selection> s_njet_two;
    std::unique_ptr<Selection> s_njet_six;
    std::unique_ptr<Selection> s_bjet_none;
    std::unique_ptr<Selection> s_bjet_one;
    std::unique_ptr<Selection> s_btight_none;
    std::unique_ptr<Selection> s_btight_one;
}; 


InvPreselection::InvPreselection(Context & ctx){
  year = extract_year(ctx);
  is_mc = ctx.get("dataset_type") == "MC";
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
  s_njet_six.reset(new NJetSelection(6));
  s_njet_two.reset(new NJetSelection(2));
  s_metFilters.reset(new METFilters(ctx));
  s_bjet_none.reset(new NJetSelection(-1,0,bmedium));
  s_bjet_one.reset(new NJetSelection(1,1,bmedium));
  s_btight_none.reset(new NJetSelection(-1,0,btight));
  s_btight_one.reset(new NJetSelection(1,1,btight));

  // Handles
  handle_event_weight = ctx.declare_event_output<double>("event_weight");
  handle_tight_b = ctx.declare_event_output<int>("tight_b");

  // Common
  common_modules.reset(new CommonModules());
  common_modules->disable_jec();
  common_modules->disable_jersmear();
  common_modules->disable_jetpfidfilter();
  common_modules->set_jet_id(jet_id);
  common_modules->switch_jetlepcleaner(true);
  common_modules->switch_jetPtSorter();
  common_modules->init(ctx);

  // Histograms
  h_baseline.reset(new AtoZHHists(ctx, "CutFlow_Baseline"));
  h_six_jets.reset(new AtoZHHists(ctx, "CutFlow_SixJets"));
  h_no_leptons.reset(new AtoZHHists(ctx, "CutFlow_LeptonVeto"));
  h_met_cut.reset(new AtoZHHists(ctx, "CutFlow_MissingPT"));
  h_delta_cut.reset(new AtoZHHists(ctx, "CutFlow_deltaphi"));

  //h_bjet_none.reset(new AtoZHHists(ctx, "NoBJets"));
  //h_bjet_one.reset(new AtoZHHists(ctx, "OneBJets"));
  h_bjet_two.reset(new AtoZHHists(ctx, "TwoBJets"));  
  h_btag_eff.reset(new BTagMCEfficiencyHists(ctx, "2_BTagMCEff", bmedium));
  h_unc_norm.reset(new NormalisationHists(ctx, "UncNorms"));

  rochester.reset(new RochesterCorrections(ctx));
  clnr_jetpuid.reset(new JetCleaner(ctx, JetPUID(JetPUID::WP_TIGHT)));
}

bool InvPreselection::process(Event & event) {

  if (is_mc) h_unc_norm->fill(event);

  // Apply rochester corrections before ID'ing muons
  rochester->process(event);
  
  //Lepton Cleaning
  clean_collection(*event.muons, event, muonId_loose);
  clean_collection(*event.electrons, event, electronId_loose);
  
  // Jet PU ID, to be applied on JEC corrected jets
  // JECs in common_modules, hence here PU ID
  bool pileup_cleaned = clnr_jetpuid->process(event);
  if (!pileup_cleaned) {return false;}

  // Apply Jet corrections
  bool passes_common = common_modules->process(event);
  if (!passes_common) { return false; }

  h_baseline->fill(event);

  //Lepton Veto  
  if ((*event.muons).size()!=0){return false;}
  if ((*event.electrons).size()!=0){return false;}
  h_no_leptons->fill(event);

  // Cut on missing transvers momentum
  double met = event.met->pt();
  if ( met<50 ) { return false; }
  h_met_cut->fill(event);  

  // MET isolation cut  
  bool has_two_jets = s_njet_two->passes(event);
  if (!has_two_jets) { return false; }
  double phi_lead = event.jets->at(0).phi();
  double phi_sub = event.jets->at(1).phi();
  double phi_miss = event.met->phi();
  double diff_lead = abs(phi_lead - phi_miss) < M_PI ? abs(phi_lead - phi_miss): 2*M_PI-abs(phi_lead - phi_miss);
  double diff_sub = abs(phi_sub - phi_miss) < M_PI ? abs(phi_sub - phi_miss): 2*M_PI-abs(phi_sub - phi_miss);
  if (diff_lead>2 && diff_sub<1) {return false;}
  if (diff_sub>2 && diff_lead<1) {return false;}
  h_delta_cut->fill(event);

  // Jet Selection
  bool has_six_jets = s_njet_six->passes(event);
  if (!has_six_jets) { return false; }
  h_six_jets->fill(event);

  // Histogram for BTagging efficiencye
  if ( !event.isRealData ) { h_btag_eff->fill(event); }

  //BTagging
  bool has_nob = s_bjet_none->passes(event);
  bool has_oneb = s_bjet_one->passes(event);
  if (has_nob || has_oneb) {return false;}
  h_bjet_two->fill(event);

  //Number of tight btags
  bool no_tight_b = s_btight_none->passes(event);
  bool one_tight_b = s_btight_one->passes(event);
  int tight_b = 2;
  if (no_tight_b) tight_b = 0;
  if (one_tight_b) tight_b = 1;
  event.set(handle_tight_b,tight_b);

  //double b_score = event.jets->at(0).btag_DeepFlavour_b()

  event.set(handle_event_weight, event.weight);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(InvPreselection)

