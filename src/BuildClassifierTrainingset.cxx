#include <math.h>

#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/Particle.h"

#include "TH1F.h"

#include "UHH2/common/include/TTbarGen.h"


using namespace std;
using namespace uhh2;


class BuildDatasetHists: public Hists {

  public:
    BuildDatasetHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
      book<TH1F>("NMatchedCombinationsPerEvent", "N_{Combs,matched}", 10, 0, 10);
      book<TH1F>("DuplicateMatches", "duplicate_matches", 2, 0, 2);
    };  
    virtual void fill(const Event & event) override { (void)event; }; 
    virtual void fill(int n_matched, int duplicate_matches) {
      hist("NMatchedCombinationsPerEvent")->Fill(n_matched);
      hist("DuplicateMatches")->Fill(duplicate_matches);
    };  
    virtual ~BuildDatasetHists() {}; 

};

class RecoMatchingDuplicatesHist: public Hists {

  public:
    RecoMatchingDuplicatesHist(Context & ctx, const string & dirname): Hists(ctx, dirname){
      book<TH1F>("RecoDuplicate", "N_{Combs,matched}", 2, 0, 2);
    };  
    virtual void fill(const Event & event) override { (void)event; }; 
    virtual void fill(int duplicate_matched) {
      hist("RecoDuplicate")->Fill(duplicate_matched);
    };  
    virtual ~RecoMatchingDuplicatesHist() {}; 
};




class BuildClassifierTrainingset: public AnalysisModule {
  public:
    explicit BuildClassifierTrainingset(Context & ctx);
    virtual bool process(Event & event) override;

    void GetJetCombinations(std::vector<Jet> comb, std::vector<Jet> jets, int k, std::vector<std::vector<Jet>> & combinations);
    bool BuildTrainingset(uhh2::Event & event);
    bool JetInVector(const Jet & j, vector<const Jet*> & combination);
    bool JetInVector(const Jet & j, vector<const GenJet*> & combination);
    double ChiSquared(Jet b1, Jet q11, Jet q12, Jet b2, Jet q21, Jet q22);

  private:
    // Common
    std::unique_ptr<CommonModules> common;

    // Histograms
    std::unique_ptr<BuildDatasetHists> hist_dataset;
    std::unique_ptr<RecoMatchingDuplicatesHist> hist_reco_matching_duplicates;

    //Selections
    std::unique_ptr<Selection> s_njet;

    // Handles
    uhh2::Event::Handle<std::vector<int>> handle_combination_labels;
    uhh2::Event::Handle<std::vector<double>> handle_chi2_score;

    uhh2::Event::Handle<std::vector<double>> handle_jet_one_pt;
    uhh2::Event::Handle<std::vector<double>> handle_jet_two_pt;
    uhh2::Event::Handle<std::vector<double>> handle_jet_thr_pt;
    uhh2::Event::Handle<std::vector<double>> handle_jet_fou_pt;
    uhh2::Event::Handle<std::vector<double>> handle_jet_fiv_pt;
    uhh2::Event::Handle<std::vector<double>> handle_jet_six_pt;

    uhh2::Event::Handle<std::vector<double>> handle_jet_one_eta;
    uhh2::Event::Handle<std::vector<double>> handle_jet_two_eta;
    uhh2::Event::Handle<std::vector<double>> handle_jet_thr_eta;
    uhh2::Event::Handle<std::vector<double>> handle_jet_fou_eta;
    uhh2::Event::Handle<std::vector<double>> handle_jet_fiv_eta;
    uhh2::Event::Handle<std::vector<double>> handle_jet_six_eta;

    uhh2::Event::Handle<std::vector<double>> handle_jet_one_phi;
    uhh2::Event::Handle<std::vector<double>> handle_jet_two_phi;
    uhh2::Event::Handle<std::vector<double>> handle_jet_thr_phi;
    uhh2::Event::Handle<std::vector<double>> handle_jet_fou_phi;
    uhh2::Event::Handle<std::vector<double>> handle_jet_fiv_phi;
    uhh2::Event::Handle<std::vector<double>> handle_jet_six_phi;
    
    uhh2::Event::Handle<std::vector<double>> handle_jet_one_e;
    uhh2::Event::Handle<std::vector<double>> handle_jet_two_e;
    uhh2::Event::Handle<std::vector<double>> handle_jet_thr_e;
    uhh2::Event::Handle<std::vector<double>> handle_jet_fou_e;
    uhh2::Event::Handle<std::vector<double>> handle_jet_fiv_e;
    uhh2::Event::Handle<std::vector<double>> handle_jet_six_e;
};


BuildClassifierTrainingset::BuildClassifierTrainingset(Context & ctx) {
  // Cleaners and Common
  const ElectronId electronId = AndId<Electron>(PtEtaSCCut(20, 2.4), ElectronID_Fall17_tight_noIso);
  const MuonId muonId = AndId<Muon>(PtEtaCut(20, 2.4), MuonID(Muon::CutBasedIdLoose), MuonIso(0.15));
  const JetId jetId = AndId<Jet>(PtEtaCut(20, 2.4), JetPFID(JetPFID::WP_TIGHT_LEPVETO));
  common.reset(new CommonModules());
  common->set_muon_id(muonId);
  common->set_jet_id(jetId);
  common->set_electron_id(electronId);
  common->switch_jetPtSorter();
  common->init(ctx);

  // Handles
  handle_combination_labels = ctx.declare_event_output<vector<int>>("jet_combination_labels");
  handle_chi2_score = ctx.declare_event_output<vector<double>>("jet_chi2_score");

  handle_jet_one_pt = ctx.declare_event_output<vector<double>>("jet_one_pt");
  handle_jet_two_pt = ctx.declare_event_output<vector<double>>("jet_two_pt");
  handle_jet_thr_pt = ctx.declare_event_output<vector<double>>("jet_thr_pt");
  handle_jet_fou_pt = ctx.declare_event_output<vector<double>>("jet_fou_pt");
  handle_jet_fiv_pt = ctx.declare_event_output<vector<double>>("jet_fiv_pt");
  handle_jet_six_pt = ctx.declare_event_output<vector<double>>("jet_six_pt");

  handle_jet_one_eta = ctx.declare_event_output<vector<double>>("jet_one_eta");
  handle_jet_two_eta = ctx.declare_event_output<vector<double>>("jet_two_eta");
  handle_jet_thr_eta = ctx.declare_event_output<vector<double>>("jet_thr_eta");
  handle_jet_fou_eta = ctx.declare_event_output<vector<double>>("jet_fou_eta");
  handle_jet_fiv_eta = ctx.declare_event_output<vector<double>>("jet_fiv_eta");
  handle_jet_six_eta = ctx.declare_event_output<vector<double>>("jet_six_eta");

  handle_jet_one_phi = ctx.declare_event_output<vector<double>>("jet_one_phi");
  handle_jet_two_phi = ctx.declare_event_output<vector<double>>("jet_two_phi");
  handle_jet_thr_phi = ctx.declare_event_output<vector<double>>("jet_thr_phi");
  handle_jet_fou_phi = ctx.declare_event_output<vector<double>>("jet_fou_phi");
  handle_jet_fiv_phi = ctx.declare_event_output<vector<double>>("jet_fiv_phi");
  handle_jet_six_phi = ctx.declare_event_output<vector<double>>("jet_six_phi");

  handle_jet_one_e = ctx.declare_event_output<vector<double>>("jet_one_e");
  handle_jet_two_e = ctx.declare_event_output<vector<double>>("jet_two_e");
  handle_jet_thr_e = ctx.declare_event_output<vector<double>>("jet_thr_e");
  handle_jet_fou_e = ctx.declare_event_output<vector<double>>("jet_fou_e");
  handle_jet_fiv_e = ctx.declare_event_output<vector<double>>("jet_fiv_e");
  handle_jet_six_e = ctx.declare_event_output<vector<double>>("jet_six_e");

  // Selections
  s_njet.reset(new NJetSelection(6));

  // Histograms
  hist_dataset.reset(new BuildDatasetHists(ctx, "Histograms"));
  hist_reco_matching_duplicates.reset(new RecoMatchingDuplicatesHist(ctx, "RecoMatchingDuplicates"));
}


void BuildClassifierTrainingset::GetJetCombinations(vector<Jet> comb, vector<Jet> jets, int k, vector<vector<Jet>> & combinations) {
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


bool BuildClassifierTrainingset::JetInVector(const Jet & j, vector<const Jet*> & combination){
  assert(combination.size() == 6);
  double lowest_delta_r = std::numeric_limits<double>::infinity();
  for ( auto & c_jet: combination ) {
    double dr = deltaR(j.v4(), c_jet->v4());
    if (dr < lowest_delta_r) lowest_delta_r = dr;
  }
  return (lowest_delta_r < .4);
}

bool BuildClassifierTrainingset::JetInVector(const Jet & j, vector<const GenJet*> & combination){
  assert(combination.size() == 6);
  double lowest_delta_r = std::numeric_limits<double>::infinity();
  for ( auto & c_jet: combination ) {
    double dr = deltaR(j.v4(), c_jet->v4());
    if (dr < lowest_delta_r) lowest_delta_r = dr;
  }
  return (lowest_delta_r < .4);
}


double BuildClassifierTrainingset::ChiSquared(Jet b1, Jet q11, Jet q12, Jet b2, Jet q21, Jet q22) {
  // Calibration values
  const double top_mass = 155.14;
  const double w_mass = 62.05;
  const double sigma_top = 34.99; 
  const double sigma_w = 13.996;

  // Compute masses
  double m_t1 = (b1.v4() + q11.v4() + q12.v4()).M();
  double m_w1 = (q11.v4() + q12.v4()).M();
  double m_t2 = (b2.v4() + q21.v4() + q22.v4()).M();
  double m_w2 = (q21.v4() + q22.v4()).M();

  // Compute Chi2 Score
  double chi_sq = pow((m_t1 - top_mass) / sigma_top, 2)
    + pow((m_t2 - top_mass) / sigma_top, 2)
    + pow((m_w1 - w_mass) / sigma_w, 2)
    + pow((m_w2 - w_mass) / sigma_w, 2); 
  return chi_sq;
}


bool BuildClassifierTrainingset::BuildTrainingset(Event & event) {
  // Get the six gen level partons and assert is hadronic decay
  TTbarGen ttbar_gen = TTbarGen(*event.genparticles);
  if (!(ttbar_gen.IsTopHadronicDecay() && ttbar_gen.IsAntiTopHadronicDecay())) {
    cout << "Decay not Hadronic" << endl;
    return false;
  }

  // Match Gen level partons to reco jets
  const GenJet* b1 = closestParticle(ttbar_gen.bTop(), *event.genjets);
  const GenJet* q11 = closestParticle(ttbar_gen.Wdecay1(), *event.genjets);
  const GenJet* q12 = closestParticle(ttbar_gen.Wdecay2(), *event.genjets);
  const GenJet* b2 = closestParticle(ttbar_gen.bAntitop(), *event.genjets);
  const GenJet* q21 = closestParticle(ttbar_gen.WMinusdecay1(), *event.genjets);
  const GenJet* q22 = closestParticle(ttbar_gen.WMinusdecay2(), *event.genjets);
  vector<const GenJet*> matched_genjets {b1, q11, q12, b2, q21, q22};

  // if the same genjet was matched to two genparticles, discard event
  for ( unsigned i = 0; i < matched_genjets.size(); i++ ){
    for ( unsigned j = i + 1; j < matched_genjets.size(); j++ ){
      if ( matched_genjets.at(i) == matched_genjets.at(j) ) {
        // hist_dataset->fill(-1, 1);
        hist_reco_matching_duplicates->fill(1);
        return false;
      }
    }
  }
  hist_reco_matching_duplicates->fill(0);

  // Build Reco Jet Combinations for Event
  int jets_per_combination = 6;
  vector<vector<Jet>> jet_combinations;
  vector<Jet> comb;
  GetJetCombinations(comb, *event.jets, jets_per_combination, jet_combinations);

  // Assert there are no duplicates in any of the reco jet_combinations
  for ( auto x: jet_combinations ) {
    for ( unsigned i = 0; i < x.size(); i++ ){
      for ( unsigned j = i + 1; j < x.size(); j++ ){
        if ( x.at(i) == x.at(j) ) {
          cout << "DUPLICATE IN JET COMBINATIONS" << endl << endl;
          return false;
        }
      }
    }
  }

  // Construct Permutations of Combinations
  int permutations_per_combination = tgamma(6 + 1); 
  int no_of_permutations = permutations_per_combination * jet_combinations.size();
  vector<vector<Jet>> ttbar_jet_permutations(no_of_permutations);

  for (unsigned i = 0; i < jet_combinations.size(); i++) {
    vector<int> v {0, 1, 2, 3, 4, 5}; 
    unsigned j = 0;
    do {
      ttbar_jet_permutations.at(permutations_per_combination*i+j).push_back(jet_combinations.at(i).at(v.at(0)));
      ttbar_jet_permutations.at(permutations_per_combination*i+j).push_back(jet_combinations.at(i).at(v.at(1)));
      ttbar_jet_permutations.at(permutations_per_combination*i+j).push_back(jet_combinations.at(i).at(v.at(2)));
      ttbar_jet_permutations.at(permutations_per_combination*i+j).push_back(jet_combinations.at(i).at(v.at(3)));
      ttbar_jet_permutations.at(permutations_per_combination*i+j).push_back(jet_combinations.at(i).at(v.at(4)));
      ttbar_jet_permutations.at(permutations_per_combination*i+j).push_back(jet_combinations.at(i).at(v.at(5)));
      j++;
    } while(next_permutation(v.begin(), v.end()));
  }

  // Find labels
  vector<int> labels (ttbar_jet_permutations.size(), 0);
  // > Iterate through permutations of reco jets
  for ( unsigned k = 0; k < ttbar_jet_permutations.size(); k++ ) {
    bool all_jets_in_perm_have_gen_match = true;
    // >> Iterate through jets in permutation
    for ( unsigned j=0; j<6; j++) {
      Jet jet = ttbar_jet_permutations.at(k).at(j);
      double dr = deltaR(jet.v4(), matched_genjets.at(j)->v4());
      bool jet_matches_gen = (dr < .4);
      if (!jet_matches_gen) {
        all_jets_in_perm_have_gen_match = false;
        break;
      }
    }
    if (all_jets_in_perm_have_gen_match) { labels.at(k) = 1; }
  }

  assert(labels.size() == ttbar_jet_permutations.size());
  // Set Handles
  event.set(handle_combination_labels, labels);
  vector<double> jet_one_pt, jet_two_pt, jet_thr_pt, jet_fou_pt, jet_fiv_pt, jet_six_pt;
  vector<double> jet_one_eta, jet_two_eta, jet_thr_eta, jet_fou_eta, jet_fiv_eta, jet_six_eta;
  vector<double> jet_one_phi, jet_two_phi, jet_thr_phi, jet_fou_phi, jet_fiv_phi, jet_six_phi;
  vector<double> jet_one_e, jet_two_e, jet_thr_e, jet_fou_e, jet_fiv_e, jet_six_e;
  vector<double> chi2_scores;
  for (unsigned i=0; i<ttbar_jet_permutations.size(); i++){
    jet_one_pt.push_back(ttbar_jet_permutations.at(i).at(0).v4().Pt());
    jet_two_pt.push_back(ttbar_jet_permutations.at(i).at(1).v4().Pt());
    jet_thr_pt.push_back(ttbar_jet_permutations.at(i).at(2).v4().Pt());
    jet_fou_pt.push_back(ttbar_jet_permutations.at(i).at(3).v4().Pt());
    jet_fiv_pt.push_back(ttbar_jet_permutations.at(i).at(4).v4().Pt());
    jet_six_pt.push_back(ttbar_jet_permutations.at(i).at(5).v4().Pt());

    jet_one_eta.push_back(ttbar_jet_permutations.at(i).at(0).v4().Eta());
    jet_two_eta.push_back(ttbar_jet_permutations.at(i).at(1).v4().Eta());
    jet_thr_eta.push_back(ttbar_jet_permutations.at(i).at(2).v4().Eta());
    jet_fou_eta.push_back(ttbar_jet_permutations.at(i).at(3).v4().Eta());
    jet_fiv_eta.push_back(ttbar_jet_permutations.at(i).at(4).v4().Eta());
    jet_six_eta.push_back(ttbar_jet_permutations.at(i).at(5).v4().Eta());

    jet_one_phi.push_back(ttbar_jet_permutations.at(i).at(0).v4().Phi());
    jet_two_phi.push_back(ttbar_jet_permutations.at(i).at(1).v4().Phi());
    jet_thr_phi.push_back(ttbar_jet_permutations.at(i).at(2).v4().Phi());
    jet_fou_phi.push_back(ttbar_jet_permutations.at(i).at(3).v4().Phi());
    jet_fiv_phi.push_back(ttbar_jet_permutations.at(i).at(4).v4().Phi());
    jet_six_phi.push_back(ttbar_jet_permutations.at(i).at(5).v4().Phi());

    jet_one_e.push_back(ttbar_jet_permutations.at(i).at(0).v4().E());
    jet_two_e.push_back(ttbar_jet_permutations.at(i).at(1).v4().E());
    jet_thr_e.push_back(ttbar_jet_permutations.at(i).at(2).v4().E());
    jet_fou_e.push_back(ttbar_jet_permutations.at(i).at(3).v4().E());
    jet_fiv_e.push_back(ttbar_jet_permutations.at(i).at(4).v4().E());
    jet_six_e.push_back(ttbar_jet_permutations.at(i).at(5).v4().E());

    // Chi2 Score
    const Jet* reco_b1 = closestParticle(*b1, ttbar_jet_permutations.at(i));
    const Jet* reco_q11 = closestParticle(*q11, ttbar_jet_permutations.at(i));
    const Jet* reco_q12 = closestParticle(*q12, ttbar_jet_permutations.at(i));
    const Jet* reco_b2 = closestParticle(*b2, ttbar_jet_permutations.at(i));
    const Jet* reco_q21 = closestParticle(*q21, ttbar_jet_permutations.at(i));
    const Jet* reco_q22 = closestParticle(*q22, ttbar_jet_permutations.at(i));
    double score = ChiSquared(*reco_b1, *reco_q11, *reco_q12, *reco_b2, *reco_q21, *reco_q22);
    chi2_scores.push_back(score);
  }
  event.set(handle_chi2_score, chi2_scores);

  event.set(handle_jet_one_pt, jet_one_pt);
  event.set(handle_jet_two_pt, jet_two_pt);
  event.set(handle_jet_thr_pt, jet_thr_pt);
  event.set(handle_jet_fou_pt, jet_fou_pt);
  event.set(handle_jet_fiv_pt, jet_fiv_pt);
  event.set(handle_jet_six_pt, jet_six_pt);

  event.set(handle_jet_one_eta, jet_one_eta);
  event.set(handle_jet_two_eta, jet_two_eta);
  event.set(handle_jet_thr_eta, jet_thr_eta);
  event.set(handle_jet_fou_eta, jet_fou_eta);
  event.set(handle_jet_fiv_eta, jet_fiv_eta);
  event.set(handle_jet_six_eta, jet_six_eta);

  event.set(handle_jet_one_phi, jet_one_phi);
  event.set(handle_jet_two_phi, jet_two_phi);
  event.set(handle_jet_thr_phi, jet_thr_phi);
  event.set(handle_jet_fou_phi, jet_fou_phi);
  event.set(handle_jet_fiv_phi, jet_fiv_phi);
  event.set(handle_jet_six_phi, jet_six_phi);

  event.set(handle_jet_one_e, jet_one_e);
  event.set(handle_jet_two_e, jet_two_e);
  event.set(handle_jet_thr_e, jet_thr_e);
  event.set(handle_jet_fou_e, jet_fou_e);
  event.set(handle_jet_fiv_e, jet_fiv_e);
  event.set(handle_jet_six_e, jet_six_e);

  int sum_of_labels = 0;
  for (int n : labels) sum_of_labels += n;
  hist_dataset->fill(sum_of_labels, 0);
  if (sum_of_labels != 1) return false;
  return true;
}


bool BuildClassifierTrainingset::process(Event & event){
  bool commonResult = common->process(event);
  if (!commonResult) return false;

  bool has_n_jets = s_njet->passes(event);
  if (!has_n_jets) return false;

  return BuildTrainingset(event);
}

UHH2_REGISTER_ANALYSIS_MODULE(BuildClassifierTrainingset)

