#ifndef _EPKE_EPKE_SOLVER_HEADER_
#define _EPKE_EPKE_SOLVER_HEADER_

#include <memory>

#include "parareal/precursor.hpp"
#include "parareal/solver_parameters.hpp"
#include "parareal/definitions.hpp"
#include "pugi/pugixml.hpp"

namespace epke {

class Solver {
public:
  template <typename T>
  using precBins    = para::precBins<T>; // binning over precursor groups
  using timeBins    = para::timeBins;    // binning over time variable
  using timeIndex   = para::timeIndex;
  using precIndex   = para::precIndex;
  using ptr         = std::shared_ptr<Solver>;
  using FineSolvers = std::vector<Solver::ptr>;

protected:
  // input member variables
  EPKEParameters::ptr params;

  // precomputed values of power, rho, and concentrations (so the simulation
  // can start at t > 0)
  EPKEOutput::ptr precomp;

  // vector of pointers to fine solvers created by this coarse solver
  FineSolvers _fine_solvers;

  // output member variables
  timeBins power;                    // power at various points in time
  timeBins rho;                      // reactivity with feedback
  precBins<timeBins> concentrations; // precursor concentrations

  // private methods
  const double computeOmega(const precIndex k, const timeIndex n) const;

  const double computeZetaHat(const precIndex k, const timeIndex n) const;

  const double computePower(const timeIndex n, const double alpha) const;

  const double computeA1(const timeIndex n) const;

  const double computeB1(const timeIndex n) const;

  const double computeDT(const timeIndex n) const;

  const double computeGamma(const timeIndex n) const;

  const bool acceptTransformation(const timeIndex n, const double alpha) const;

  virtual const double getPower(const timeIndex n) const { return power.at(n); }

  virtual const double getConcentration(const precIndex k,
					const timeIndex n) const {
    return concentrations.at(k).at(n);
  }

  virtual const timeIndex getFirstTimeIndex() const {
    return precomp->getNumTimeSteps();
  }

public:
  // Construct from pugixml nodes
  Solver(const pugi::xml_node& input_node, const pugi::xml_node& output_node);

  // Construct from parameters and precomputed objects
  Solver(const EPKEParameters::ptr parameters,
	 const EPKEOutput::ptr precomputed);

  // Construct a fine solver by interpolating parameters from coarse solver
  Solver::ptr createFineSolver(const timeBins& fine_time,
			       const timeIndex coarse_index);


  para::SolverOutput::ptr updateCoarseSolution();

  // Assemble global output from vector of fine solver output
  const EPKEOutput::ptr assembleGlobalOutput() const;

  // TODO: Make this a virtual function from the Solver base class
  EPKEOutput::ptr solve();

  EPKEParameters::ptr getParameters() const {
    return params;
  }

  EPKEOutput::ptr getPrecomputed() const {
    return precomp;
  }

  EPKEOutput::ptr getOutput() const {
    return std::make_shared<EPKEOutput>(params->getTime(),
					concentrations,
					power,
					rho);
  }

  const timeIndex getNumPrecompTimeSteps() const {
    return precomp->getNumTimeSteps();
  }
};

class Propagator : public Solver {
private:
  const double getPower(const timeIndex n) const override {
    return precomp->getPower(n);
  }

  const double getConcentration(const precIndex k,
				const timeIndex n) const override {
    return precomp->getConcentration(k,n);
  }

  const timeIndex getFirstTimeIndex() const override { return 1; }

public:
  Propagator(const pugi::xml_node& input_node,
	     const pugi::xml_node& output_node) :
    Solver(input_node, output_node) {}

  Propagator(const EPKEParameters::ptr parameters,
	     const EPKEOutput::ptr     precomputed) :
    Solver(parameters, precomputed) {}

};

} // namespace epke

#endif
