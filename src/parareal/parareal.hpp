#ifndef _PARAREAL_PARAREAL_HEADER_
#define _PARAREAL_PARAREAL_HEADER_

#include "parareal/definitions.hpp"
#include "parareal/solver_output.hpp"
#include "parareal/solver_parameters.hpp"
#include "epke/epke_solver.hpp"

// TODO: Write the solver base class and make the EPKE solver a derived class

namespace para {

class Parareal {
public:
  using timeBins  = para::timeBins;
  using timeIndex = para::timeIndex;
  using precIndex = para::precIndex;

private:
  // Max number of parareal iterations
  const precIndex _max_iterations;

  // Number of fine time steps for every coarse time step
  const timeIndex _n_fine_per_coarse;

  // Fine time grid (low fidelity) solver
  epke::Solver _solver;

public:
  Parareal(const epke::Solver&     solver,
	   const precIndex         max_iterations,
	   const timeIndex         n_fine_per_coarse)
    : _parameters(parameters),
      _output(output),
      _solver(solver),
      _max_iterations(max_iterations),
      _n_fine_per_coarse(n_fine_per_coarse) {}

  // Generate fine time
  timeBins generateFineTime(const timeIndex n);

  // Update the coarse solution with a parareal iteration
  timeBins update_coarse_solution();

  // Solve
  void solve() { _solver.solve(); }

};
} // namespace para

#endif
