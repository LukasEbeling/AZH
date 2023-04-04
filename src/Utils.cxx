#include "UHH2/core/include/Event.h"
#include "UHH2/AZH/include/Utils.h"


using namespace std;
using namespace uhh2;


PDGUtils::PDGUtils() {}


bool PDGUtils::IsMuon(const GenParticle* p) {
  int id = p->pdgId();
  return abs(id) == 13;
}


bool PDGUtils::IsMuon(GenParticle& p) {
  int id = p.pdgId();
  return abs(id) == 13;
}


bool PDGUtils::IsElectron(GenParticle& p) {
  int id = p.pdgId();
  return abs(id) == 11;
}


bool PDGUtils::IsW(const GenParticle* p) {
  if (!p) return false;
  int id = p->pdgId();
  return abs(id) == 24;
}


bool PDGUtils::IsZ(const GenParticle* p) {
  if (!p) return false;
  int id = p->pdgId();
  return abs(id) == 23;
}


namespace pdgIdUtils {

  bool is_Muon(const GenParticle & gp) {
    int id = gp.pdgId();
    return abs(id) == 13;
  }

  bool is_H(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 35;
  }

  bool is_A(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 36;
  }

  bool is_top(const GenParticle & gp){
    int id = gp.pdgId();
    return id == 6;
  }

  bool is_antitop(const GenParticle & gp){
    int id = gp.pdgId();
    return id == -6;
  }

  bool is_Z(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 23;
  }

  bool is_W(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 24;
  }

  bool is_WPlus(const GenParticle & gp){
    int id = gp.pdgId();
    return id == 24;
  }

  bool is_WMinus(const GenParticle & gp){
    int id = gp.pdgId();
    return id == -24;
  }

  bool is_b_bbar(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 5;
  }

  bool is_bbar(const GenParticle & gp){
    int id = gp.pdgId();
    return id == -5;
  }

  bool is_b(const GenParticle & gp){
    int id = gp.pdgId();
    return id == 5;
  }

  bool is_light_quark(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 1 || id == 2 || id == 3 || id == 4;
  }

  bool is_emu(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 11 || id == 13;
  }

  bool is_charged_lepton(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 11 || id == 13 || id == 15;
  }

  bool is_neutrino(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 12 || id == 14 || id == 16;
  }

}


// Taken from https://github.com/MatthiesC/LegacyTopTagging/blob/master/src/Utils.cxx
// Included muons into filtering, as done in central UHH2 module and as per recommendation.
HEMSelection::HEMSelection(Context& ctx){
  fYear = extract_year(ctx);
}

bool HEMSelection::passes(const Event & event){
  if(fYear == Year::isUL18 && ((event.isRealData && event.run >= fRunNumber) || !event.isRealData)){
    vector<Particle> relevant_objects;
    relevant_objects.insert(relevant_objects.end(), event.jets->begin(), event.jets->end());
    relevant_objects.insert(relevant_objects.end(), event.electrons->begin(), event.electrons->end());
    relevant_objects.insert(relevant_objects.end(), event.muons->begin(), event.muons->end());
    for(const Particle & obj : relevant_objects) {
      if(obj.v4().eta() > fEtaRange.first && obj.v4().eta() < fEtaRange.second && obj.v4().phi() > fPhiRange.first && obj.v4().phi() < fPhiRange.second) {
        return true;
      }
    }
  }
  return false;
}