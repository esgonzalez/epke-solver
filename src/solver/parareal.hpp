#ifndef _PARAREAL_PARAREAL_HEADER_
#define _PARAREAL_PARAREAL_HEADER_

#include "solver_parameters.hpp"

// TODO: Write the solver base class and make the EPKE solver a derived class
//#include "solver.hpp"

class Parareal {
private:
  using timeBins = std::vector<double>; // binning over time variable
  
  // Max number of parareal iterations
  int _K;

  // Number of fine time steps for every coarse time step
  int _n_fine_per_coarse;

  // Time dependent solution on the coarse time grid
  timeBins _coarse_time;
  timeBins _coarse_solution;

  // Time dependent solver parameters on the coarse time grid. These will be
  // interpolated for the fine solver
  SolverParameters _coarse_parameters;

  // Fine time grid (low fidelity) solver
  //Solver _fine_solver;
  
public:
  Parareal(timeBins coarse_time,
	   timeBins coarse_solution,
	   int K,
	   int n_fine_per_coarse)
    : _coarse_time(coarse_time),
      _coarse_solution(coarse_solution),
      _K(K),
      _n_fine_per_coarse(n_fine_per_coarse) {}

  // Update the coarse solution with a parareal iteration
  timeBins update_coarse_solution();
  
};

#endif
