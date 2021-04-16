#ifndef _PARAREAL_SOLVER_PARAMETERS_HEADER_
#define _PARAREAL_SOLVER_PARAMETERS_HEADER_

#include <vector>

#include "parareal/precursor.hpp"

class SolverParameters {
protected:
  template <typename T>
  using precBins = std::vector<T>;      // binning over precursor groups
  using timeBins = std::vector<double>; // binning over time variable
  using timeIndex = uint32_t;
  using precIndex = uint8_t;
};

class EPKEParameters : SolverParameters {
private:
  // coarse time mesh parameters
  const timeBins _time;     // vector of time points
  const timeBins _rho_imp;  // the imposed reactivity (without feedback)
  const timeBins _gen_time; // Mean neutron generation time (Lambda)
  const timeBins _pow_norm; // Power normalization factor
  const timeBins _beta_eff; // total delayed neutron fraction
  const timeBins _lambda_h; // Linear heat conduction constant
  const precBins<Precursor::ptr> _precursors;
  
public:
  EPKEParameters(const timeBins& time,
		 const timeBins& rho_imp,
		 const timeBins& gen_time,
		 const timeBins& pow_norm,		 
		 const timeBins& beta_eff,
		 const timeBins& lambda_h,
		 const precBins<Precursor::ptr> precursors)
    : SolverParameters(),
      _time(time),
      _rho_imp(rho_imp),
      _gen_time(gen_time),
      _pow_norm(pow_norm),
      _beta_eff(beta_eff),
      _lambda_h(lambda_h),      
      _precursors(precursors) {}

  // Getters
  const precIndex getNumPrecursors() const { return _precursors.size(); }
  const timeIndex getNumTimeSteps()  const { return _time.size(); } 
  
  const timeBins& getTime()    const { return _time;     }
  const timeBins& getGenTime() const { return _gen_time; }
  const timeBins& getPowNorm() const { return _pow_norm; }
  const timeBins& getBetaEff() const { return _beta_eff; }
  const timeBins& getLambdaH() const { return _lambda_h; }

  const double getTime(timeIndex n)   const { return _time.at(n); }
  const double getRhoImp(timeIndex n) const { return _rho_imp.at(n); }
  const double getGenTime(timeIndex n) const { return _gen_time.at(n); }
  
  const precBins<Precursor::ptr>& getPrecursors() const { return _precursors; }
  Precursor::ptr getPrecursor(precIndex k) const { return _precursors.at(k); }

  const timeBins& getDelayedFraction(precIndex k) const {
    return _precursors.at(k)->delayedFraction();
  }

  const double getDelayedFraction(precIndex k, timeIndex n) const {
    return _precursors.at(k)->delayedFraction(n);
  }

  const timeBins& getDecayConstant(precIndex k) const {
    return _precursors.at(k)->decayConstant();
  }

  const double getDecayConstant(precIndex k, timeIndex n) const {
    return _precursors.at(k)->decayConstant(n);
  }

  
  // Interpolate parameters for the fine time mesh
  const EPKEParameters interpolate(const timeBins& fine_time) const;
  
};

#endif
