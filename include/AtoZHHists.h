#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"

#include "UHH2/AZH/include/Utils.h"


/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */

class AtoZHHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    AtoZHHists(uhh2::Context & ctx, const std::string & dirname);

    // Methods
    virtual void fill(const uhh2::Event & event) override;
    virtual ~AtoZHHists();
};
