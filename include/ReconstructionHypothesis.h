#pragma once

#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/TopJet.h"
#include <map>

/**
 * Abstraction of the different ttbar decays that can be puzzled together from
 * the jet constituents of the event.
 */
class ReconstructionHypothesis {
  private:
    // Members
    LorentzVector t_b_v4;
    LorentzVector t_q1_v4;
    LorentzVector t_q2_v4;
    LorentzVector tbar_b_v4;
    LorentzVector tbar_q1_v4;
    LorentzVector tbar_q2_v4;

  public:
    // LorentzVector Getters
    LorentzVector get_t_b_v4() const { return t_b_v4; };
    LorentzVector get_t_q1_v4() const { return t_q1_v4; };
    LorentzVector get_t_q2_v4() const { return t_q2_v4; };
    LorentzVector get_tbar_b_v4() const { return tbar_b_v4; };
    LorentzVector get_tbar_q1_v4() const { return tbar_q1_v4; };
    LorentzVector get_tbar_q2_v4() const { return tbar_q2_v4; };
    LorentzVector get_t_v4() const { return t_b_v4 + t_q1_v4 + t_q2_v4; };
    LorentzVector get_t_w_v4() const { return t_q1_v4 + t_q2_v4; };
    LorentzVector get_tbar_w_v4() const { return tbar_q1_v4 + tbar_q2_v4; };
    LorentzVector get_tbar_v4() const { return tbar_b_v4 + tbar_q1_v4 + tbar_q2_v4; };
    LorentzVector get_H_v4() const { return get_t_v4() + get_tbar_v4(); };

    // Setter
    void set_t_b_v4(LorentzVector v4) { t_b_v4=v4; };
    void set_t_q1_v4(LorentzVector v4) { t_q1_v4=v4; };
    void set_t_q2_v4(LorentzVector v4) { t_q2_v4=v4; };
    void set_tbar_b_v4(LorentzVector v4) { tbar_b_v4=v4; };
    void set_tbar_q1_v4(LorentzVector v4) { tbar_q1_v4=v4; };
    void set_tbar_q2_v4(LorentzVector v4) { tbar_q2_v4=v4; };

    // Complex Getter
    double get_t_m() { return get_t_v4().M(); };
    double get_t_w_m() const { return get_t_w_v4().M(); };
    double get_tbar_w_m() { return get_tbar_w_v4().M(); };
    double get_tbar_m() { return get_tbar_v4().M(); };
    double get_H_m() { return get_H_v4().M(); };
    LorentzVector get_H_v4() { return t_b_v4 + t_q1_v4 + t_q2_v4 + tbar_b_v4 + tbar_q1_v4 + tbar_q2_v4; };

    // ~ReconstructionHypothesis();
};

