#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/AZH/include/ReconstructionHypothesis.h"


class ReconstructionHypothesisBuilder {
  public:
    ReconstructionHypothesisBuilder(uhh2::Context & ctx);

    uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> handle_ttbar_reco_hypotheses_all;
    uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> handle_ttbar_reco_hypotheses_fixedb;

    void BuildRecoHypotheses(uhh2::Event & event);

  private:
    void GetJetCombinations(std::vector<Jet> comb, std::vector<Jet> jets, int k, std::vector<std::vector<Jet>> & combinations);
    void SplitEventJetsIntoBAndLightJets(std::vector<Jet> & event_jets, std::vector<Jet> & b_jets, std::vector<Jet> & light_jets);
    void BuildAllRecoHypotheses(uhh2::Event & event);
    void BuildFixedBRecoHypotheses(uhh2::Event & event);
};

