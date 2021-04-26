#ifndef _PARAREAL_PARAREAL_HEADER_
#define _PARAREAL_PARAREAL_HEADER_

#include <fstream>

#include "parareal/definitions.hpp"
#include "parareal/solver_output.hpp"
#include "parareal/solver_parameters.hpp"
#include "pugi/pugixml.hpp"
#include "epke/epke_solver.hpp"

// TODO: Write the solver base class and make the EPKE solver a derived class

namespace para {

class Parareal {
public:
  using timeBins  = para::timeBins;
  using timeIndex = para::timeIndex;
  using precIndex = para::precIndex;

private:
  // Fine time grid (low fidelity) solver
  epke::Solver _solver;

  // Output path
  std::string _outpath;

  // Max number of parareal iterations
  const precIndex _max_iterations;

  // Number of fine time steps for every coarse time step
  const timeIndex _n_fine_per_coarse;

public:
  Parareal(const pugi::xml_node& parareal_node);

  Parareal(const epke::Solver&     solver,
	   const std::string       outpath,
	   const precIndex         max_iterations,
	   const timeIndex         n_fine_per_coarse)
    : _solver(solver),
      _outpath(outpath),
      _max_iterations(max_iterations),
      _n_fine_per_coarse(n_fine_per_coarse) {}

  // Get outpath
  std::string getOutpath() const { return _outpath; }

  // Generate fine time
  timeBins generateFineTime(const timeIndex n);

  // Update the coarse solution with a parareal iteration
  timeBins update_coarse_solution();

  // Solve
  void solve();

  // Build xml output doc
  void buildXMLDoc(pugi::xml_document& doc) const;

};
} // namespace para

#endif
