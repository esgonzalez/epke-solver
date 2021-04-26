#ifndef _PARAREAL_SOLVER_OUTPUT_HEADER_
#define _PARAREAL_SOLVER_OUTPUT_HEADER_

#include "parareal/definitions.hpp"
#include "pugi/pugixml.hpp"

namespace para {

class SolverOutput {
public:
  template <typename T>
  using precBins  = para::precBins<T>;
  using timeBins  = para::timeBins;
  using timeIndex = para::timeIndex;
  using precIndex = para::precIndex;
  using ptr       = std::shared_ptr<SolverOutput>;
};
} // namespace para

namespace epke {
class EPKEOutput : public para::SolverOutput {
private:
  // time-dependent parameters
  const timeBins _power;                    // reactor power
  const timeBins _rho;                      // reactivity with feedback
  const precBins<timeBins> _concentrations; // precursor concentrations;

public:
  // Construct from pugixml node
  EPKEOutput(const pugi::xml_node& output_node);

  // Construct from data vectors
  EPKEOutput(const timeBins& power,
	     const timeBins& rho,
	     const precBins<timeBins>& concentrations)
    : SolverOutput(),
      _power(power),
      _rho(rho),
      _concentrations(concentrations) {}

  const timeIndex getNumTimeSteps()  const { return _power.size();          }
  const precIndex getNumPrecursors() const { return _concentrations.size(); }

  const double getPower(const timeIndex n) const { return _power.at(n); }
  const double getRho(const timeIndex n)   const { return _rho.at(n);   }
  const double getConcentration(const precIndex k, const timeIndex n) const {
    return _concentrations.at(k).at(n);
  }

  const EPKEOutput createPrecomputed(const timeIndex n) const;

};
} // namespace epke

#endif
