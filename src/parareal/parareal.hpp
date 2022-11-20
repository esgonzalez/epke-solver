#ifndef _PARAREAL_PARAREAL_HEADER_
#define _PARAREAL_PARAREAL_HEADER_

#include <memory>

#include "parareal/definitions.hpp"
#include "parareal/solver_output.hpp"
#include "parareal/solver_parameters.hpp"

namespace para {

  template <typename Coarse, typename Fine>
  class Parareal {
  public:
    template <typename T>
    using precBins  = para::precBins<T>; // binning over precursor groups
    using timeBins  = para::timeBins;    // binning over time variable
    using timeIndex = para::timeIndex;
    using precIndex = para::precIndex;
    using paraIndex = para::paraIndex;
    using Output    = typename Coarse::Output;
    using Params    = typename Coarse::Params;

  private:
    // Output path
    std::string _outpath;

    // Global fine output
    typename Output::ptr _global_output;

    // Max number of parareal iterations
    const precIndex _max_iterations;

    // Coarse solver
    typename Coarse::ptr _coarse_solver;

    // Fine solver
    typename Fine::ptr _fine_solver;

    // Number of fine time steps per coarse time step
    const timeIndex _n_fine_per_coarse;

    // Result from coarse solver at parareal iteration k-1
    typename Output::ptr _old_coarse;

    // Result from coarse solver at parareal iteration k
    typename Output::ptr _new_coarse;

    // Result from the fine solvers
    typename Output::ptr _fine_coarsened;

  public:
    Parareal(typename Coarse::ptr       coarse_solver,
	     typename Fine::ptr         fine_solver,
	     typename Fine::Output::ptr global_output,
	     const    timeIndex         n_fine_per_coarse,
	     const    paraIndex         max_iterations,
	     const    std::string       outpath)
      : _coarse_solver(coarse_solver),
	_fine_solver(fine_solver),
	_n_fine_per_coarse(n_fine_per_coarse),
	_max_iterations(max_iterations),
	_outpath(outpath),
	_global_output(global_output),
	_old_coarse(std::make_shared<Output>(*_coarse_solver->getSolution())),
	_new_coarse(std::make_shared<Output>(*_coarse_solver->getSolution())) {}

    // Generate the vector of fine time steps given the coarse time steps
    timeBins generateFineTime(const timeIndex n) const;

    // Solve the coarse propagator
    void runCoarseSolver();

    // Solve the fine solver at coarse time index n
    void runFineSolver(const timeIndex n);

    // Update the solution with a parareal iteration
    void update(const paraIndex k);

     // Assemble global output from vector of fine solver output
    void assembleGlobalOutput();

    // Get outpath
    std::string getOutpath() const { return _outpath; }

    // Solve
    void solve();

    // Build xml output doc
    void writeToXML(pugi::xml_document& doc) const;
  };

} // namespace para

#include "parareal/parareal.t.hpp"

#endif
