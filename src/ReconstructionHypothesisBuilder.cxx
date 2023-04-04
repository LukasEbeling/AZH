#include <math.h>

#include "UHH2/AZH/include/ReconstructionHypothesisBuilder.h"
#include "UHH2/AZH/include/ReconstructionHypothesis.h"


using namespace std;
using namespace uhh2;


ReconstructionHypothesisBuilder::ReconstructionHypothesisBuilder(Context & ctx) {
  handle_ttbar_reco_hypotheses_all = ctx.get_handle<vector<ReconstructionHypothesis>>("ttbar_reco_hypothese_all");
  handle_ttbar_reco_hypotheses_fixedb = ctx.get_handle<vector<ReconstructionHypothesis>>("ttbar_reco_hypothese_fixedb");
}


void ReconstructionHypothesisBuilder::GetJetCombinations(vector<Jet> comb, vector<Jet> jets, int k, vector<vector<Jet>> & combinations) {
  if (k == 1) {
    for (unsigned i = 0; i < jets.size(); i++) {
      vector<Jet> _comb = comb;
      _comb.push_back(jets.at(i));
      combinations.push_back(_comb);
    }   
  }
  else {
    for (unsigned i = 0; i < jets.size(); i++) {
      vector<Jet> _comb = comb;
      _comb.push_back(jets.at(i));
      vector<Jet> _remaining_jets = jets;
      _remaining_jets.erase(_remaining_jets.begin(), _remaining_jets.begin()+i+1);
      GetJetCombinations(_comb, _remaining_jets, k-1, combinations);
    }   
  }
}


void ReconstructionHypothesisBuilder::BuildAllRecoHypotheses(Event & event) {
  vector<Jet> high_pt_reco_jets = *event.jets;
  if (high_pt_reco_jets.size() > 9) { high_pt_reco_jets.resize(9); }

  int jets_per_combination = 6;
  int permutations_per_combination = tgamma(jets_per_combination + 1);
  assert (permutations_per_combination == 720);

  vector<vector<Jet>> combinations;
  vector<Jet> comb;
  GetJetCombinations(comb, *event.jets, jets_per_combination, combinations);

  int no_of_hypotheses = permutations_per_combination * combinations.size();
  vector<ReconstructionHypothesis> all_reco_hypotheses(no_of_hypotheses);

  for (unsigned i = 0; i < combinations.size(); i++) {
    vector<int> v {0, 1, 2, 3, 4, 5};
    unsigned j = 0;
    do {
      all_reco_hypotheses.at(permutations_per_combination*i+j).set_t_b_v4(combinations.at(i).at(v.at(0)).v4());
      all_reco_hypotheses.at(permutations_per_combination*i+j).set_t_q1_v4(combinations.at(i).at(v.at(1)).v4());
      all_reco_hypotheses.at(permutations_per_combination*i+j).set_t_q2_v4(combinations.at(i).at(v.at(2)).v4());
      all_reco_hypotheses.at(permutations_per_combination*i+j).set_tbar_b_v4(combinations.at(i).at(v.at(3)).v4());
      all_reco_hypotheses.at(permutations_per_combination*i+j).set_tbar_q1_v4(combinations.at(i).at(v.at(4)).v4());
      all_reco_hypotheses.at(permutations_per_combination*i+j).set_tbar_q2_v4(combinations.at(i).at(v.at(5)).v4());
      j++;
    } while(next_permutation(v.begin(), v.end()));
  }

  event.set(handle_ttbar_reco_hypotheses_all, all_reco_hypotheses);
}


void ReconstructionHypothesisBuilder::SplitEventJetsIntoBAndLightJets(vector<Jet> & event_jets, vector<Jet> & b_jets, vector<Jet> & light_jets) {
  unsigned idx_b1 = 0;
  unsigned idx_b2 = 1;
  double score_b1 = event_jets.size();
  double score_b2 = event_jets.size();
  for ( unsigned i=2; i<event_jets.size(); i++ ) {
    double score = event_jets.at(i).btag_DeepJet();
    if ( score_b1 < score ) {
      score_b1 = score;
      idx_b1 = i;
    } else if ( score_b2 < score ) {
      score_b2 = score;
      idx_b2 = i;
    }
  }

  for ( unsigned i = 0; i<event_jets.size(); i++ ) {
    if ( i == idx_b1 || i == idx_b2 ) {
      b_jets.push_back(event_jets.at(i));
    } else { light_jets.push_back(event_jets.at(i)); }
  }
}


void ReconstructionHypothesisBuilder::BuildFixedBRecoHypotheses(Event & event) {
  vector<Jet> high_pt_reco_jets = *event.jets;
  if (high_pt_reco_jets.size() > 9) { high_pt_reco_jets.resize(9); }

  int jets_per_combination = 4;
  int permutations_per_combination = 2 * tgamma(jets_per_combination + 1);
  assert (permutations_per_combination == 48);

  // find highest btagged jets
  vector<Jet> b_jets, light_jets;
  SplitEventJetsIntoBAndLightJets(high_pt_reco_jets, b_jets, light_jets);

  vector<vector<Jet>> combinations;
  vector<Jet> comb;
  GetJetCombinations(comb, light_jets, jets_per_combination, combinations);

  int no_of_hypotheses = permutations_per_combination * combinations.size();
  vector<ReconstructionHypothesis> fixedb_reco_hypotheses(no_of_hypotheses);

  for (unsigned i = 0; i < combinations.size(); i++) {
    vector<int> v1 {0, 1};
    vector<int> v2 {0, 1, 2, 3};
    unsigned j = 0;
    do {
      do {
        fixedb_reco_hypotheses.at(permutations_per_combination*i+j).set_t_b_v4(b_jets.at(v1.at(0)).v4());
        fixedb_reco_hypotheses.at(permutations_per_combination*i+j).set_t_q1_v4(combinations.at(i).at(v2.at(0)).v4());
        fixedb_reco_hypotheses.at(permutations_per_combination*i+j).set_t_q2_v4(combinations.at(i).at(v2.at(1)).v4());
        fixedb_reco_hypotheses.at(permutations_per_combination*i+j).set_tbar_b_v4(b_jets.at(v1.at(1)).v4());
        fixedb_reco_hypotheses.at(permutations_per_combination*i+j).set_tbar_q1_v4(combinations.at(i).at(v2.at(2)).v4());
        fixedb_reco_hypotheses.at(permutations_per_combination*i+j).set_tbar_q2_v4(combinations.at(i).at(v2.at(3)).v4());
        j++;
      } while(next_permutation(v2.begin(), v2.end()));
    } while(next_permutation(v1.begin(), v1.end()));
  }

  event.set(handle_ttbar_reco_hypotheses_fixedb, fixedb_reco_hypotheses);
}


void ReconstructionHypothesisBuilder::BuildRecoHypotheses(Event & event) {
  BuildAllRecoHypotheses(event);
  BuildFixedBRecoHypotheses(event);
}

