#include "UHH2/core/include/Event.h"

#include "UHH2/AZH/include/JetCleaner.h"


using namespace std;
using namespace uhh2;


// https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
JetPUID::JetPUID(const wp & working_point): fWP(working_point) {}

bool JetPUID::operator()(const Jet & jet, const Event & event) const {
  const string year = event.year;
  const double eta = fabs(jet.v4().eta());
  const double pt = jet.v4().pt();
  // Apply PU ID only to jets with 10 < pT < 50 GeV
  if(eta > 5.0 || pt < 10 || pt > 50) return true;

  double x(-2.); // ranges between -1 and +1 (BDT discriminator)
  switch(fWP) {
    case WP_LOOSE :
      if(year == "UL16preVFP" || year == "UL16postVFP") {
        if(eta < 2.5) {
          if(pt < 20) x = -0.95;
          else if(pt < 30) x = -0.90;
          else if(pt < 40) x = -0.71;
          else x = -0.42;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.70;
          else if(pt < 30) x = -0.57;
          else if(pt < 40) x = -0.36;
          else x = -0.09;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.52;
          else if(pt < 30) x = -0.43;
          else if(pt < 40) x = -0.29;
          else x = -0.14;
        }
        else {
          if(pt < 20) x = -0.49;
          else if(pt < 30) x = -0.42;
          else if(pt < 40) x = -0.23;
          else x = -0.02;
        }
      }
      else if(year == "UL17" || year == "UL18") {
        if(eta < 2.5) {
          if(pt < 20) x = -0.95;
          else if(pt < 30) x = -0.88;
          else if(pt < 40) x = -0.63;
          else x = -0.19;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.72;
          else if(pt < 30) x = -0.55;
          else if(pt < 40) x = -0.18;
          else x = 0.22;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.68;
          else if(pt < 30) x = -0.60;
          else if(pt < 40) x = -0.43;
          else x = -0.13;
        }
        else {
          if(pt < 20) x = -0.47;
          else if(pt < 30) x = -0.43;
          else if(pt < 40) x = -0.24;
          else x = -0.03;
        }
      }
      else throw runtime_error((string)"JetPUID::operator()(): Year '"+year+"' not implemented");
      break;
    case WP_MEDIUM :
      if(year == "UL16preVFP" || year == "UL16postVFP") {
        if(eta < 2.5) {
          if(pt < 20) x = 0.20;
          else if(pt < 30) x = 0.62;
          else if(pt < 40) x = 0.86;
          else x = 0.93;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.56;
          else if(pt < 30) x = -0.39;
          else if(pt < 40) x = -0.10;
          else x = 0.19;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.43;
          else if(pt < 30) x = -0.32;
          else if(pt < 40) x = -0.15;
          else x = 0.04;
        }
        else {
          if(pt < 20) x = -0.38;
          else if(pt < 30) x = -0.29;
          else if(pt < 40) x = -0.08;
          else x = 0.12;
        }
      }
      else if(year == "UL17" || year == "UL18") {
        if(eta < 2.5) {
          if(pt < 20) x = 0.26;
          else if(pt < 30) x = 0.68;
          else if(pt < 40) x = 0.90;
          else x = 0.96;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.33;
          else if(pt < 30) x = -0.04;
          else if(pt < 40) x = 0.36;
          else x = 0.61;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.54;
          else if(pt < 30) x = -0.43;
          else if(pt < 40) x = -0.16;
          else x = 0.14;
        }
        else {
          if(pt < 20) x = -0.37;
          else if(pt < 30) x = -0.30;
          else if(pt < 40) x = -0.09;
          else x = 0.12;
        }
      }
      else throw runtime_error((string)"JetPUID::operator()(): Year '"+year+"' not implemented");
      break;
    case WP_TIGHT :
      if(year == "UL16preVFP" || year == "UL16postVFP") {
        if(eta < 2.5) {
          if(pt < 20) x = 0.71;
          else if(pt < 30) x = 0.87;
          else if(pt < 40) x = 0.94;
          else x = 0.97;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.32;
          else if(pt < 30) x = -0.08;
          else if(pt < 40) x = 0.24;
          else x = 0.48;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.30;
          else if(pt < 30) x = -0.16;
          else if(pt < 40) x = 0.05;
          else x = 0.26;
        }
        else {
          if(pt < 20) x = -0.22;
          else if(pt < 30) x = -0.12;
          else if(pt < 40) x = 0.10;
          else x = 0.29;
        }
      }
      else if(year == "UL17" || year == "UL18") {
        if(eta < 2.5) {
          if(pt < 20) x = 0.77;
          else if(pt < 30) x = 0.90;
          else if(pt < 40) x = 0.96;
          else x = 0.98;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = 0.38;
          else if(pt < 30) x = 0.60;
          else if(pt < 40) x = 0.82;
          else x = 0.92;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.31;
          else if(pt < 30) x = -0.12;
          else if(pt < 40) x = 0.20;
          else x = 0.47;
        }
        else {
          if(pt < 20) x = -0.21;
          else if(pt < 30) x = -0.13;
          else if(pt < 40) x = 0.09;
          else x = 0.29;
        }
      }
      else throw runtime_error((string)"JetPUID::operator()(): Year '"+year+"' not implemented");
      break;
    default :
      throw runtime_error((string)"JetPUID::operator()(): Unknown working point");
  }

  return jet.pileupID() > x;
}