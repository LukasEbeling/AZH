#pragma once
#include <iostream>
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/Utils.h"


class PDGUtils{
  public:
    PDGUtils();
    bool IsMuon(const GenParticle* p);
    bool IsMuon(GenParticle& p);
    bool IsElectron(GenParticle& p);
    bool IsW(const GenParticle* p);
    bool IsZ(const GenParticle* p);
    bool IsTopDaughter(const GenParticle* p);
};


// map between the channels and the PDG values of the
// corresponding product of lepton particle numbers
// e = 11, mu = 13
enum class Channel {
  diMuon=13*13,
  diElectron=11*11,
  ElectronMuon=11*13,
  Invalid=999,
};


enum class Region {
  SR = 0,
  CR_1L = 1,
  CR_lowdelta = 2,
  CR_0B = 3,
  CR_0B_2L = 4,
  CR_lowmet = 5,
  CR_1L_anymet = 6,
  Invalid=999,
};

namespace pdgIdUtils {
  bool is_Muon(const GenParticle & gp);
  bool is_H(const GenParticle & gp);
  bool is_A(const GenParticle & gp);
  bool is_top(const GenParticle & gp);
  bool is_antitop(const GenParticle & gp);
  bool is_Z(const GenParticle & gp);
  bool is_W(const GenParticle & gp);
  bool is_WPlus(const GenParticle & gp);
  bool is_WMinus(const GenParticle & gp);
  bool is_b_bbar(const GenParticle & gp);
  bool is_bbar(const GenParticle & gp);
  bool is_b(const GenParticle & gp);
  bool is_light_quark(const GenParticle & gp);
  bool is_emu(const GenParticle & gp);
  bool is_charged_lepton(const GenParticle & gp);
  bool is_neutrino(const GenParticle & gp);
}


// HEM Issue: remove events from HE problematic sectors
// https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
class HEMSelection: public uhh2::Selection {
public:
  HEMSelection(uhh2::Context & ctx);
  virtual bool passes(const uhh2::Event & event) override;
  double GetAffectedLumiFraction() const { return fAffectedLumiFraction; };
private:
  Year fYear;
  const int fRunNumber = 319077;
  const std::pair<double, double> fEtaRange = {-3.2, -1.3};
  const std::pair<double, double> fPhiRange = {-1.57, -0.87};
  const double fAffectedLumiFraction = 0.64844705699; // (Run 319077 (17.370008/pb) + Run C + Run D) / all 2018
};
