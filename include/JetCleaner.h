#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/Utils.h"

// https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
class JetPUID {
 public:
  enum wp {WP_LOOSE, WP_MEDIUM, WP_TIGHT};
  explicit JetPUID(const wp & working_point);
  bool operator()(const Jet & jet, const uhh2::Event & event) const;
 private:
  const wp fWP;
};
