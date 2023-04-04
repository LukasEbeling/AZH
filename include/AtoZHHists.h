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


class AtoZHRecoHists: public uhh2::Hists {
  public:
    AtoZHRecoHists(uhh2::Context & ctx, const std::string & dirname);

    // Handles
    uhh2::Event::Handle<double> handle_z_mass_reco;
    uhh2::Event::Handle<double> handle_z_pt_reco;
    uhh2::Event::Handle<double> handle_z_mass_gen_matched;
    uhh2::Event::Handle<double> handle_z_pt_gen_matched;
    uhh2::Event::Handle<double> handle_h_mass_reco_chi_sq;
    uhh2::Event::Handle<double> handle_t_mass_reco_chi_sq;
    uhh2::Event::Handle<double> handle_tbar_mass_reco_chi_sq;
    uhh2::Event::Handle<double> handle_h_mass_reco_six_jets;
    uhh2::Event::Handle<double> handle_h_mass_gen_six_jets;

    uhh2::Event::Handle<double> handle_a_mass;
    uhh2::Event::Handle<double> handle_a_minus_h_mass;

    // Methods
    virtual void fill(const uhh2::Event & event) override;
    ~AtoZHRecoHists();
};


class AtoZHZMassWindowScanHists: public uhh2::Hists {
  public:
    AtoZHZMassWindowScanHists(uhh2::Context & ctx, const std::string & dirname);

    // Handles
    uhh2::Event::Handle<double> handle_z_mass_reco;

    // Methods
    virtual void fill(const uhh2::Event & event) override;
    ~AtoZHZMassWindowScanHists();
};


class StandardModelBRHists: public uhh2::Hists {
  public:
    StandardModelBRHists(uhh2::Context & ctx, const std::string & dirname);

    // Methods
    virtual void fill(const uhh2::Event & event) override;
    ~StandardModelBRHists();
};

class TriggerHists: public uhh2::Hists {
  public:
    TriggerHists(uhh2::Context & ctx, const std::string & dirname);
    uhh2::Event::Handle<int> handle_channel;
    uhh2::Event::Handle<std::vector<Muon>> handle_muons_tight;
    uhh2::Event::Handle<std::vector<Electron>> handle_electrons_tight;

    // Methods
    virtual void fill(const uhh2::Event & event) override;
    ~TriggerHists();
};

