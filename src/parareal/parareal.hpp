#ifndef _PARAREAL_PARAREAL_HEADER_
#define _PARAREAL_PARAREAL_HEADER_

#include "parareal/solver_parameters.hpp"

// TODO: Write the solver base class and make the EPKE solver a derived class
//#include "solver.hpp"

class Parareal {
public:
  using TimeBins  = std::vector<double>; // binning over time variable
  using TimeIndex = uint32_t;

private:
  // Max number of parareal iterations
  int _K;

  // Number of fine time steps for every coarse time step
  int _n_fine_per_coarse;

  // Time dependent solution on the coarse time grid
  TimeBins _coarse_time;
  TimeBins _coarse_solution;

  // Time dependent solver parameters on the coarse time grid. These will be
  // interpolated for the fine solver
  SolverParameters _coarse_parameters;

  // Fine time grid (low fidelity) solver
  //Solver _fine_solver;
  
public:
  Parareal(TimeBins coarse_time,
	   TimeBins coarse_solution,
	   int K,
	   int n_fine_per_coarse)
    : _coarse_time(coarse_time),
      _coarse_solution(coarse_solution),
      _K(K),
      _n_fine_per_coarse(n_fine_per_coarse) {}

  // Generate fine time
  TimeBins generateFineTime(TimeIndex n);
  
  // Update the coarse solution with a parareal iteration
  TimeBins update_coarse_solution();
  
};

#endif
