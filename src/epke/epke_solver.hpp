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

private:
  // input member variables
  const EPKEParameters params;

  // precomputed values of power, rho, and concentrations (so the simulation
  // can start at t > 0)
  const EPKEOutput precomp;

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

public:
  // Construct from pugixml nodes
  Solver(const pugi::xml_node& input_node, const pugi::xml_node& output_node);

  // Construct from parameters and precomputed objects
  Solver(const EPKEParameters& parameters, const EPKEOutput& precomputed);

  // Construct a fine solver by interpolating parameters from coarse solver
  Solver::ptr createFineSolver(const timeBins& fine_time,
			       const timeIndex coarse_index);

  // Assemble global output from vector of fine solver output
  const para::SolverOutput::ptr assembleGlobalOutput() const;

  // TODO: Make this a virtual function from the Solver base class
  para::SolverOutput::ptr solve();

  para::SolverParameters::ptr getParameters() const {
    return std::make_shared<EPKEParameters>(params);
  }

  para::SolverOutput::ptr getPrecomputed() const {
    return std::make_shared<EPKEOutput>(precomp);
  }

  const timeIndex getNumPrecompTimeSteps() const {
    return precomp.getNumTimeSteps();
  }
};
} // namespace epke

#endif
