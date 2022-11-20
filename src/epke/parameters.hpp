#ifndef _EPKE_PARAMETERS_HEADER_
#define _EPKE_PARAMETERS_HEADER_

#include <memory>

#include "parareal/solver_parameters.hpp"
#include "epke/precursor.hpp"

namespace epke {

class EPKEParameters : public para::SolverParameters {
public:
  template <typename T>
  using precBins  = para::precBins<T>;
  using precIndex = para::precIndex;
  using Base      = SolverParameters;
  using ptr       = std::shared_ptr<EPKEParameters>;

private:
  // time-dependent parameters
  const timeBins _rho_imp;  // imposed reactivity (without feedback)
  const timeBins _gen_time; // mean neutron generation time (Lambda)
  const timeBins _pow_norm; // power normalization factor
  const timeBins _beta_eff; // total delayed neutron fraction
  const timeBins _lambda_h; // linear heat conduction constant

  // Pointers to precursor objects
  const precBins<Precursor::ptr> _precursors;

  // coefficient for finite differencing scheme
  // theta = 0 -> fully explicit, theta = 1 -> fully implicit
  const double _theta;

  // linear thermal feedback coefficient
  const double _gamma_d;

  // eta = 0 -> first order heat conduction for total power
  // eta = 1 -> first order heat conduction for power increment
  const double _eta;

public:
  EPKEParameters(const timeBins& time,
		 const precBins<Precursor::ptr>& precursors,
		 const timeBins& rho_imp,
		 const timeBins& gen_time,
		 const timeBins& pow_norm,
		 const timeBins& beta_eff,
		 const timeBins& lambda_h,
		 const double theta,
		 const double gamma_d,
		 const double eta)
    : SolverParameters(time),
      _precursors(precursors),
      _rho_imp(rho_imp),
      _gen_time(gen_time),
      _pow_norm(pow_norm),
      _beta_eff(beta_eff),
      _lambda_h(lambda_h),
      _theta(theta),
      _gamma_d(gamma_d),
      _eta(eta) {}

  // Getters
  const double getTheta()                    const { return _theta;           }
  const double getGammaD()                   const { return _gamma_d;         }
  const double getEta()                      const { return _eta;             }
  const double getRhoImp(const timeIndex n)  const { return _rho_imp.at(n);   }
  const double getGenTime(const timeIndex n) const { return _gen_time.at(n);  }
  const double getPowNorm(const timeIndex n) const { return _pow_norm.at(n);  }
  const double getBetaEff(const timeIndex n) const { return _beta_eff.at(n);  }
  const double getLambdaH(const timeIndex n) const { return _lambda_h.at(n);  }

  const precIndex getNumPrecursors() const override {
    return _precursors.size();
  }

  const double getDelayedFraction(precIndex k, timeIndex n) const {
    return _precursors.at(k)->delayedFraction(n);
  }

  const double getDecayConstant(precIndex k, timeIndex n) const {
    return _precursors.at(k)->decayConstant(n);
  }

  // Interpolate parameters for the fine time mesh
  virtual Base::ptr interpolateImpl(const timeBins& fine_time) const override;

  void writeToXML(pugi::xml_document& doc) const override;
};

} // namespace epke

#endif
