#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/PSWeights.h"

#include "UHH2/AZH/include/NormalisationTools.h"
#include "UHH2/AZH/include/ScaleFactors.h"
#include "UHH2/AZH/include/Utils.h"//for region
#include "UHH2/AZH/include/HiggsReco.h"
#include "UHH2/AZH/include/AtoZHHists.h"


using namespace std;
using namespace uhh2;

class InvReconstruction : public AnalysisModule{
  public:
  explicit InvReconstruction(Context& ctx);
  virtual bool process(Event& event) override;

  private:
  //Methods
  bool assign_region(Event& event);

  //Handles
  Event::Handle<int> handle_region;
  Event::Handle<double> handle_sum_pt;
  Event::Handle<double> handle_weight;

  unique_ptr<HiggsReconstructor> higgs_reconstructor;
  //unique_ptr<AnalysisModule> pdf_weights;
  unique_ptr<HEMSelection> sel_hem;
  unique_ptr<MCLumiWeight> sf_lumi;
  unique_ptr<MCPileupReweight> sf_pileup;
  unique_ptr<AnalysisModule> sf_l1prefiring;
  unique_ptr<AnalysisModule> sf_vjets;
  unique_ptr<AnalysisModule> sf_mtop;
  unique_ptr<PSWeights> ps_weights;
  unique_ptr<MCScaleVariation> sf_QCDScaleVar;

  std::unique_ptr<Hists> h_kinematics;
};

InvReconstruction::InvReconstruction(Context& ctx){
  handle_region = ctx.declare_event_output<int>("region");
  handle_sum_pt = ctx.declare_event_output<double>("HT");
  handle_weight = ctx.get_handle<double>("event_weight");

  higgs_reconstructor.reset(new HiggsReconstructor(ctx));
  //pdf_weights.reset(new PDFWeightHandleProducer(ctx));
  sel_hem.reset(new HEMSelection(ctx));
  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_pileup.reset(new MCPileupReweight(ctx));
  sf_l1prefiring.reset(new L1PrefiringWeight(ctx));
  sf_vjets.reset(new VJetsReweighting(ctx));
  sf_mtop.reset(new TopPtReweighting(ctx, string2bool(ctx.get("apply_TopPtReweighting"))));
  sf_QCDScaleVar.reset(new MCScaleVariation(ctx));

  h_kinematics.reset(new RecoHists(ctx, "Kinematics"));
}

bool InvReconstruction::process(Event& event){
  event.weight = event.get(handle_weight);//Wozu? 
  double weight = event.weight;

  //pdf_weights->process(event);
  if(sel_hem->passes(event)) {
    if(event.isRealData) return false;
    else event.weight *= (1. - sel_hem->GetAffectedLumiFraction());
  }
  sf_lumi->process(event);
  bool PUProfileCheck = sf_pileup->process(event);
  sf_l1prefiring->process(event);
  sf_vjets->process(event);
  sf_mtop->process(event);
  //ps_weights->process(event); not working yet
  bool LHEWeightsCheck = sf_QCDScaleVar->process(event);
  if (!event.isRealData && (!PUProfileCheck && !LHEWeightsCheck)) { throw std::invalid_argument("Failed checks on PU and LHE weights."); }
  //To Do: Btag efficiency

  bool valid_region = assign_region(event);
  if(!valid_region) return false;

  event.weight = weight;
  higgs_reconstructor->process(event);

  double HT = 0;
  for(Jet jet: *event.jets){HT += jet.pt();}
  event.set(handle_sum_pt,HT);

  event.set(handle_weight, event.weight);

  return true;
}

bool InvReconstruction::assign_region(Event& event){
  Region region = Region::Invalid;

  //bool has_six_jets = s_njet_six->passes(event);
  //bool has_no_b = passes_minmax(*event.jets, 0, 0, event, btag_id_one);
  //bool has_one_b = passes_minmax(*event.jets, 1, 1, event, btag_id_one);
  //bool has_two_b = passes_minmax(*event.jets, 2, -1, event, btag_id_two);
  //bool met_window = (0 <= event.met->pt())&&(event.met->pt()<=300);

  //if(has_six_jets && has_two_b && met_window) region = Region::SignalRegion;
  region = Region::SignalRegion;
  event.set(handle_region, (int) region);
  h_kinematics->fill(event);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(InvReconstruction)