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
    double GetTransMassA();
    double GetTransMassH();
    double GetMassH();
    double GetTransMassDiff();
    double GetMassDiff();

private:
    unique_ptr<ReconstructionHypothesisBuilder> reconstruction_hypothesis_builder;
    Event::Handle<vector<ReconstructionHypothesis>> handle_ttbar_reco_hypotheses_fixedb;
    
    //Constants
    const double top_mass = 164.79;
    const double w_mass = 80.00;
    const double sigma_top = 27.60; 
    const double sigma_w = 16.86;
    const double z_mass = 91.19;

    //Higgs
    double mt_A = 0;
    double mt_H = 0;
    double m_H = 0;

    //Methods 
    double ChiSquared(ReconstructionHypothesis hypothesis);
    void TransMassAtGen(Event & event);
    LorentzVector ProcessFiveJets(Event & event);
    LorentzVector ProcessMoreJets(Event & event);
    LorentzVector Projection(LorentzVector v);
    ReconstructionHypothesis BestHypothesis(vector<ReconstructionHypothesis> hypotheses);
};

