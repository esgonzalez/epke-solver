#ifndef _PARAREAL_SOLVER_PARAMETERS_HEADER_
#define _PARAREAL_SOLVER_PARAMETERS_HEADER_

#include <vector>

#include "parareal/precursor.hpp"
#include "parareal/definitions.hpp"
#include "pugi/pugixml.hpp"

namespace para {

class SolverParameters {
public:
  template <typename T>
  using precBins  = para::precBins<T>;
  using timeBins  = para::timeBins;
  using timeIndex = para::timeIndex;
  using precIndex = para::precIndex;

protected:
  const timeBins _time;                       // time points
  const precBins<Precursor::ptr> _precursors; // pointers to precursor objects

private:
  static const precBins<Precursor::ptr>
  buildPrecursors(const pugi::xml_node& precursors_node,
		  const timeIndex n_steps) {
    precBins<Precursor::ptr> precursors;
    for (const auto& prec_node : precursors_node) {
      precursors.push_back(std::make_shared<Precursor>(prec_node, n_steps));
    }
    return precursors;
  }

public:
  // Construct from pugixml node
  SolverParameters(const pugi::xml_node& solver_node);

  // Construct from data vectors
  SolverParameters(const timeBins& time,
		   const precBins<Precursor::ptr>& precursors)
    : _time(time), _precursors(precursors) {}

  // Getters
  const precIndex getNumPrecursors()   const { return _precursors.size(); }
  const timeIndex getNumTimeSteps()    const { return _time.size();       }
  const double    getTime(timeIndex n) const { return _time.at(n);        }

  const double getDelayedFraction(precIndex k, timeIndex n) const {
    return _precursors.at(k)->delayedFraction(n);
  }

  const double getDecayConstant(precIndex k, timeIndex n) const {
    return _precursors.at(k)->decayConstant(n);
  }

};
} // namespace para

// TODO: Break up this class into its own hpp and cpp files and move to epke
//       directory
class EPKEParameters : public para::SolverParameters {
private:
  // time-dependent parameters
  const timeBins _rho_imp;  // imposed reactivity (without feedback)
  const timeBins _gen_time; // mean neutron generation time (Lambda)
  const timeBins _pow_norm; // power normalization factor
  const timeBins _beta_eff; // total delayed neutron fraction
  const timeBins _lambda_h; // linear heat conduction constant

  // precomputed output values to allow simulate to start at t > 0
  const timeBins           _p_history;      // power history
  const precBins<timeBins> _zeta_histories; // concentration histories

  // coefficient for finite differencing scheme
  // theta = 0 -> fully explicit, theta = 1 -> fully implicit
  const double _theta;

  // linear thermal feedback coefficient
  const double _gamma_d;

  // eta = 0 -> first order heat conduction for total power
  // eta = 1 -> first order heat conduction for power increment
  const double _eta;

public:
  // Construct from pugi xml node
  EPKEParameters(const pugi::xml_node& epke_node);

  // Construct from data vectors
  EPKEParameters(const timeBins& time,
		 const precBins<Precursor::ptr>& precursors,
		 const timeBins& rho_imp,
		 const timeBins& gen_time,
		 const timeBins& pow_norm,
		 const timeBins& beta_eff,
		 const timeBins& lambda_h,
		 const timeBins& p_history,
		 const precBins<timeBins>& zeta_histories,
		 const double theta,
		 const double gamma_d,
		 const double eta)
    : SolverParameters(time, precursors),
      _rho_imp(rho_imp),
      _gen_time(gen_time),
      _pow_norm(pow_norm),
      _beta_eff(beta_eff),
      _lambda_h(lambda_h),
      _p_history(p_history),
      _zeta_histories(zeta_histories),
      _theta(theta),
      _gamma_d(gamma_d),
      _eta(eta) {}

  // Getters
  const timeIndex getNumPrecomputedTimeSteps() const {
    return _p_history.size();
  }

  const double getTheta()  const { return _theta;   }
  const double getGammaD() const { return _gamma_d; }
  const double getEta()    const { return _eta;     }

  const double getRhoImp(timeIndex n)   const { return _rho_imp.at(n);   }
  const double getGenTime(timeIndex n)  const { return _gen_time.at(n);  }
  const double getPowNorm(timeIndex n)  const { return _pow_norm.at(n);  }
  const double getBetaEff(timeIndex n)  const { return _beta_eff.at(n);  }
  const double getLambdaH(timeIndex n)  const { return _lambda_h.at(n);  }

  // These getters are for grabbing precomputed values from previous simulations
  const double getPHistory(timeIndex n) const { return _p_history.at(n); }

  const double getZetaHistory(precIndex k, timeIndex n) const {
    return _zeta_histories.at(k).at(n);
  }

  // TODO: The power and concentration histories for the fine mesh will have to
  //       come from the results of the coarse mesh
  // Interpolate parameters for the fine time mesh
  const EPKEParameters interpolate(const timeBins& fine_time) const;

};

#endif
