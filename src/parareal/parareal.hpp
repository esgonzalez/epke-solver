#ifndef _PARAREAL_PARAREAL_HEADER_
#define _PARAREAL_PARAREAL_HEADER_

#include "parareal/definitions.hpp"
#include "parareal/solver_output.hpp"
#include "parareal/solver_parameters.hpp"

// TODO: Write the solver base class and make the EPKE solver a derived class
//#include "solver.hpp"

namespace para {

class Parareal {
public:
  using timeBins  = para::timeBins;
  using timeIndex = para::timeIndex;
  using precIndex = para::precIndex;

private:
  // Max number of parareal iterations
  const precIndex _K;

  // Number of fine time steps for every coarse time step
  const timeIndex _n_fine_per_coarse;

  // Time dependent solution on the coarse time grid
  const SolverOutput _coarse_output;

  // Time dependent solver parameters on the coarse time grid. These will be
  // interpolated for the fine solver
  const SolverParameters _coarse_parameters;

  // Fine time grid (low fidelity) solver
  //Solver _fine_solver;

public:
  Parareal(const SolverParameters& coarse_parameters,
	   const SolverOutput&     coarse_output,
	   const precIndex         K,
	   const timeIndex         n_fine_per_coarse)
    : _coarse_parameters(coarse_parameters),
      _coarse_output(coarse_output),
      _K(K),
      _n_fine_per_coarse(n_fine_per_coarse) {}

  // Generate fine time
  timeBins generateFineTime(const timeIndex n);

  // Update the coarse solution with a parareal iteration
  timeBins update_coarse_solution();

};
} // namespace para

#endif
