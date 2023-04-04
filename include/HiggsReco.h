#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"

#include "UHH2/AZH/include/ReconstructionHypothesis.h"
#include "UHH2/AZH/include/ReconstructionHypothesisBuilder.h"

using namespace uhh2;
using namespace std;

class HiggsReconstructor: uhh2::AnalysisModule {
public:
    explicit HiggsReconstructor(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & event);

    Event::Handle<double> handle_H_m;
    Event::Handle<double> handle_met;
    Event::Handle<double> handle_A_mt;
    Event::Handle<double> handle_H_mt;
    Event::Handle<double> handle_W_m;
    Event::Handle<double> handle_gen_A_mt;

private:
    unique_ptr<ReconstructionHypothesisBuilder> reconstruction_hypothesis_builder;
    Event::Handle<vector<ReconstructionHypothesis>> handle_ttbar_reco_hypotheses_fixedb;
    const double top_mass = 164.79;
    const double w_mass = 80.00;
    const double sigma_top = 27.60; 
    const double sigma_w = 16.86;
    const double z_mass = 91.19;
    double ChiSquared(ReconstructionHypothesis hypothesis);
    double TransMassAtGen(Event & event);
    double WMass(ReconstructionHypothesis hypothesis);
    double BestWMass(vector<ReconstructionHypothesis> hypotheses);
    LorentzVector ProcessFiveJets(Event & event);
    LorentzVector ProcessMoreJets(Event & event);
    LorentzVector Projection(LorentzVector v);
    ReconstructionHypothesis BestHypothesis(vector<ReconstructionHypothesis> hypotheses);
};

