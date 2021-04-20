#ifndef _PARAREAL_SOLVER_OUTPUT_HEADER_
#define _PARAREAL_SOLVER_OUTPUT_HEADER_

#include "parareal/definitions.hpp"

namespace para {

class SolverOutput {
protected:
  template <typename T>
  using precBins  = para::precBins<T>;
  using timeBins  = para::timeBins;
  using timeIndex = para::timeIndex;
  using precIndex = para::precIndex;
};
} // namespace para

namespace epke {
class EPKEOutput : para::SolverOutput {
private:
  // time-dependent parameters
  const timeBins _power;                    // reactor power
  const timeBins _rho;                      // reactivity with feedback
  const precBins<timeBins> _concentrations; // precursor concentrations;

public:
  EPKEOutput(const timeBins& power,
	     const timeBins& rho,
	     const precBins<timeBins>& concentrations)
    : SolverOutput(),
      _power(power),
      _rho(rho),
      _concentrations(concentrations) {}

  // TODO: Create a function, getHistory(n), to return an EPKEOutput object
  //       truncated after n so we don't have to call getPower(n),
  //       getConcentration(k,n) every time we want to access precomputed data
  const double getPower(const timeIndex n) const { return _power.at(n); }
  const double getConcentration(const precIndex k, const timeIndex n) const {
    return _concentrations.at(k).at(n);
  }

};
} // namespace epke

#endif
