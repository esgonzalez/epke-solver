#ifndef _PARAREAL_SOLVER_OUTPUT_HEADER_
#define _PARAREAL_SOLVER_OUTPUT_HEADER_

#include "parareal/definitions.hpp"

namespace para {

class SolverOutput {
protected:
  template <typename T>
  using precBins = para::precBins<T>;
  using timeBins = para::timeBins;
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
};
} // namespace epke

#endif
