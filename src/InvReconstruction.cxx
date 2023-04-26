#include <iostream>
#include <memory>
#include <cmath>

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/PSWeights.h"

//#include "UHH2/AZH/include/NormalisationTools.h"
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
  unique_ptr<HiggsReconstructor> higgs_reconstructor;
  
  //Methods
  bool assign_region(Event& event);
  double delta_phi(Event& event);

  //Handles
  Event::Handle<int> handle_region;
  Event::Handle<double> handle_sum_pt;
  Event::Handle<double> handle_weight;
  Event::Handle<double> handle_delta;

  //Weights missing in preselection
  std::unique_ptr<PSWeights> ps_weights;
  std::unique_ptr<AnalysisModule> pdf_weights;

  //Histogram Pointers
  std::unique_ptr<Hists> h_preselection;
  std::unique_ptr<Hists> h_met_180;
  std::unique_ptr<Hists> h_delta_phi;
};

InvReconstruction::InvReconstruction(Context& ctx){
  handle_region = ctx.declare_event_output<int>("region");
  handle_sum_pt = ctx.declare_event_output<double>("HT");
  handle_weight = ctx.get_handle<double>("event_weight");
  handle_delta = ctx.declare_event_output<double>("delta_phi");

  higgs_reconstructor.reset(new HiggsReconstructor(ctx));

  //pdf_weights.reset(new PDFWeightHandleProducer(ctx)); not working yet
  ps_weights.reset(new PSWeights(ctx));

  h_preselection.reset(new PreHists(ctx, "CutFlow_Preselection"));
  h_met_180.reset(new PreHists(ctx, "CutFlow_MET>180"));
  h_delta_phi.reset(new PreHists(ctx, "CutFlow_deltaphi"));
}

bool InvReconstruction::process(Event& event){
  event.weight = event.get(handle_weight);
  ///pdf_weights->process(event);
  ps_weights->process(event);

  h_preselection->fill(event);


  if(event.met->pt()<170) return false;
  h_met_180->fill(event);

  double delta_min = delta_phi(event);
  if(delta_min < 0.5) return false;
  event.set(handle_delta,delta_min);
  h_delta_phi->fill(event);

  bool valid_region = assign_region(event);
  if(!valid_region) return false;

  higgs_reconstructor->process(event);

  double HT = 0;
  for(Jet jet: *event.jets){HT += jet.pt();}
  event.set(handle_sum_pt,HT);

  event.set(handle_weight, event.weight);

  return true;
}

double InvReconstruction::delta_phi(Event& event){
  double phi_m = event.met->phi();
  double min_diff = M_PI;

  for (Jet jet: *event.jets){
    if(jet.pt()<30) continue;
    double phi_j = jet.phi();
    double delta = abs(phi_m - phi_j);
    if(delta > M_PI){delta = 2*M_PI - delta;}
    if(delta < min_diff) min_diff = delta;
  }

  return min_diff;
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

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(InvReconstruction)