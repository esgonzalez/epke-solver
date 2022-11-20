#ifndef _PARAREAL_SOLVER_HEADER_
#define _PARAREAL_SOLVER_HEADER_

#include "parareal/definitions.hpp"
#include "parareal/solver_parameters.hpp"
#include "parareal/solver_output.hpp"

namespace para {

  class Solver {
  public:
    template <typename T>
    using precBins  = para::precBins<T>; // binning over precursor groups
    using timeBins  = para::timeBins;    // binning over time variable
    using timeIndex = para::timeIndex;
    using precIndex = para::precIndex;
    using Params    = para::SolverParameters;
    using Output    = para::SolverOutput;

  protected:
    // Solver parameters
    Params::ptr _params;

    // Solver solution
    Output::ptr _solution;

  public:
    Solver(Params::ptr parameters, Output::ptr solution)
      : _params(parameters), _solution(solution) {}

    void baseReset(Output::ptr solution) {
      _solution = solution;
    }

    const double getTime(const timeIndex n) const {
      return _params->getTime(n);
    }

    const timeBins& getTime() const {
      return _params->getTime();
    }

    const timeIndex getNumTimeSteps() const {
      return _solution->getNumTimeSteps();
    }

    const precIndex getNumPrecursors() const {
      return _params->getNumPrecursors();
    }

    const timeIndex getStartTimeIndex() const {
      return _solution->getStartTimeIndex();
    }

    const timeIndex getStopTimeIndex() const {
      return _solution->getStopTimeIndex();
    }

    // Advance one time step
    virtual void step(const timeIndex n) = 0;

    // Solve the problem
    virtual void solve() = 0;
}; // class Solver



} // namespace para

#endif
