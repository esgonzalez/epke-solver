#ifndef _PARAREAL_SOLVER_OUTPUT_HEADER_
#define _PARAREAL_SOLVER_OUTPUT_HEADER_

#include <memory>

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

protected:
  const timeBins           _time;           // time mesh
  const precBins<timeBins> _concentrations; // precursor concentrations

public:
  // Construct from pugixml node
  SolverOutput(const pugi::xml_node & precomp_node);

  // Construct from data vectors
  SolverOutput(const timeBins& time,
	       const precBins<timeBins>& concentrations)
    : _time(time), _concentrations(concentrations) {}

  // Getters
  const timeIndex getNumTimeSteps()  const { return _time.size(); }
  const precIndex getNumPrecursors() const { return _concentrations.size(); }

  const double getConcentration(const precIndex k, const timeIndex n) const {
    return _concentrations.at(k).at(n);
  }

  // Create output object with truncated precomputed values from coarse solver
  virtual SolverOutput::ptr createPrecomputedImpl(const timeIndex n) const = 0;

  // Create output object with only information from the coarse time steps
  virtual SolverOutput::ptr coarsenImpl(const timeBins& coarse_time) const = 0;

  // Write the output to an xml document
  virtual void writeToXML(pugi::xml_document& doc) const = 0;
};

  template<typename T>
  std::shared_ptr<T> createPrecomputed(std::shared_ptr<T> precomp,
				       const timeIndex n) {
    return std::static_pointer_cast<T>(precomp->createPrecomputedImpl(n));
  }

  template<typename T>
  std::shared_ptr<T> coarsen(std::shared_ptr<T> fine_output,
			     const timeBins& coarse_time) {
    return std::static_pointer_cast<T>(fine_output->coarsenImpl(coarse_time));
  }

} // namespace para

namespace epke {

class EPKEOutput : public para::SolverOutput {
public:
  using ptr = std::shared_ptr<EPKEOutput>;

private:
  // time-dependent parameters
  const timeBins _power; // reactor power
  const timeBins _rho;   // reactivity with feedback

public:
  // Construct from pugixml node
  EPKEOutput(const pugi::xml_node& precomp_node);

  // Construct from data vectors
  EPKEOutput(const timeBins&           time,
	     const precBins<timeBins>& concentrations,
	     const timeBins&           power,
	     const timeBins&           rho)
    : SolverOutput(time, concentrations), _power(power), _rho(rho) {}

  const double getPower(const timeIndex n) const { return _power.at(n); }
  const double getRho(const timeIndex n)   const { return _rho.at(n);   }

  SolverOutput::ptr createPrecomputedImpl(const timeIndex n) const override;

  SolverOutput::ptr coarsenImpl(const timeBins& coarse_time) const override;

  void writeToXML(pugi::xml_document& doc) const override;

};

} // namespace epke

#endif
