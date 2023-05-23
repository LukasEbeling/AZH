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

#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/PSWeights.h"

#include "UHH2/AZH/include/AtoZHHists.h"
#include "UHH2/AZH/include/JetCleaner.h"
#include "UHH2/AZH/include/MetFilters.h"
#include "UHH2/AZH/include/NormalisationTools.h"
#include "UHH2/AZH/include/RochesterCorrections.h"
#include "UHH2/AZH/include/ScaleFactors.h"
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
    double angle_difference(double,double);
    bool is_mc;
    Year year;
    const JetId jet_id = AndId<Jet>(PtEtaCut(30, 2.4), JetPFID(JetPFID::WP_TIGHT_LEPVETO_CHS));
    MuonId muonId_loose;
    ElectronId electronId_loose;

    // Selection niherits from AnalysisModule as it uses passes() and process()
    std::unique_ptr<AnalysisModule> s_metFilters;

    // Handles
    uhh2::Event::Handle<double> handle_event_weight;
    uhh2::Event::Handle<double> handle_origin_weight;

    // Histograms
    std::unique_ptr<Hists> h_unc_norm;
    std::unique_ptr<Hists> h_baseline;
    std::unique_ptr<Hists> h_no_leptons;
    std::unique_ptr<Hists> h_six_jets;
    //std::unique_ptr<Hists> h_bjet_two;
    std::unique_ptr<Hists> h_met_100;

    // Histograms for BTagging efficiency measurements
    std::unique_ptr<BTagMCEfficiencyHists> h_btag_eff;

    // Rochester corrections and leptons ID+ISO collections
    std::unique_ptr<AnalysisModule> rochester;

    // Selection
    std::unique_ptr<CommonModules> common_modules;
    std::unique_ptr<AnalysisModule> clnr_jetpuid;
    std::unique_ptr<Selection> s_njet_two;
    std::unique_ptr<Selection> s_njet_six;
    std::unique_ptr<Selection> s_bjet_none;
    std::unique_ptr<Selection> s_bjet_one;
    std::unique_ptr<Selection> s_btight_none;
    std::unique_ptr<Selection> s_btight_one;

    // Event Weighting
    unique_ptr<HEMSelection> sel_hem;
    unique_ptr<MCLumiWeight> sf_lumi;
    unique_ptr<MCPileupReweight> sf_pileup;
    unique_ptr<AnalysisModule> sf_l1prefiring;
    unique_ptr<AnalysisModule> sf_vjets;
    unique_ptr<AnalysisModule> sf_mtop;
    unique_ptr<MCScaleVariation> sf_QCDScaleVar;
    unique_ptr<PSWeights> ps_weights;
    unique_ptr<AnalysisModule> pdf_weights;
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

  // Event Weighting
  sel_hem.reset(new HEMSelection(ctx));
  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_pileup.reset(new MCPileupReweight(ctx));
  sf_l1prefiring.reset(new L1PrefiringWeight(ctx));
  sf_vjets.reset(new VJetsReweighting(ctx));
  sf_mtop.reset(new TopPtReweighting(ctx, string2bool(ctx.get("apply_TopPtReweighting"))));
  sf_QCDScaleVar.reset(new MCScaleVariation(ctx));
  pdf_weights.reset(new PDFWeightHandleProducer(ctx)); 
  ps_weights.reset(new PSWeights(ctx));

  // Handles
  handle_event_weight = ctx.declare_event_output<double>("event_weight");
  handle_origin_weight = ctx.declare_event_output<double>("origin_weight");

  // Common
  common_modules.reset(new CommonModules());
  common_modules-> disable_mclumiweight();
  common_modules-> disable_mcpileupreweight();
  //common_modules-> disable_lumisel();
  //common_modules-> disable_metfilters();
  //common_modules-> disable_pvfilter();
  common_modules->disable_jec();
  common_modules->disable_jersmear();
  common_modules->disable_jetpfidfilter();
  common_modules->set_jet_id(jet_id);
  common_modules->switch_jetlepcleaner(true);
  common_modules->switch_jetPtSorter();
  common_modules->init(ctx);

  // Histograms
  h_baseline.reset(new PreHists(ctx, "CutFlow_Baseline"));
  h_six_jets.reset(new PreHists(ctx, "CutFlow_SixJets"));
  h_no_leptons.reset(new PreHists(ctx, "CutFlow_LeptonVeto"));
  h_met_100.reset(new PreHists(ctx, "CutFlow_MET>100"));
  //h_bjet_two.reset(new PreHists(ctx, "CutFlow_TwoB"));
  
  h_btag_eff.reset(new BTagMCEfficiencyHists(ctx, "2_BTagMCEff", bmedium));
  h_unc_norm.reset(new NormalisationHists(ctx, "UncNorms"));

  rochester.reset(new RochesterCorrections(ctx));
  clnr_jetpuid.reset(new JetCleaner(ctx, JetPUID(JetPUID::WP_TIGHT)));
}

bool InvPreselection::process(Event & event) {
  
  event.set(handle_origin_weight,event.weight);

  // Event Weighting
  sf_lumi->process(event);
  sf_pileup->process(event);
  sf_l1prefiring->process(event);
  sf_vjets->process(event);
  sf_mtop->process(event); 
  sf_QCDScaleVar->process(event);
  pdf_weights->process(event);
  ps_weights->process(event);
  
  if(sel_hem->passes(event)) {
    if(event.isRealData) return false;
    else event.weight *= (1. - sel_hem->GetAffectedLumiFraction());
  }  
  
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
  int leptons = (*event.electrons).size() + (*event.muons).size();
  if (leptons>2){return false;}
  h_no_leptons->fill(event);

  // Cut on missing transvers momentum
  double met = event.met->pt(); 
  if (met<100) return false;
  h_met_100->fill(event);

  // Jet Selection
  bool has_six_jets = s_njet_six->passes(event);
  if (!has_six_jets) { return false; }
  h_six_jets->fill(event);

  //BTagging
  //bool has_nob = s_bjet_none->passes(event);
  //bool has_oneb = s_bjet_one->passes(event);
  //if (has_nob || has_oneb) {return false;}
  //h_bjet_two->fill(event)

  // Histogram for BTagging efficiencye
  if ( !event.isRealData ) { h_btag_eff->fill(event); }

  //double b_score = event.jets->at(0).btag_DeepFlavour_b()

  event.set(handle_event_weight, event.weight);

  return true;
}

double InvPreselection::angle_difference(double jet_phi, double met_phi){
  double diff = abs(jet_phi-met_phi);
  if(diff > M_PI){diff = 2*M_PI - diff;}
  return diff;
}

UHH2_REGISTER_ANALYSIS_MODULE(InvPreselection)

