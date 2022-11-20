#ifndef _EPKE_SOLVER_HEADER_
#define _EPKE_SOLVER_HEADER_

#include <memory>

#include "parareal/solver.hpp"
#include "parareal/definitions.hpp"
#include "epke/parameters.hpp"
#include "epke/output.hpp"

namespace epke {

  class Solver : public para::Solver {
  public:
    using base = para::Solver;
    template <typename T>
    using precBins  = base::precBins<T>; // binning over precursor groups
    using timeBins  = base::timeBins;    // binning over time variable
    using timeIndex = base::timeIndex;
    using precIndex = base::precIndex;
    using ptr       = std::shared_ptr<epke::Solver>;
    using Output    = epke::EPKEOutput;
    using Params    = epke::EPKEParameters;

  private:
    // Input member variables
    Params::ptr _params;

    // Power, reactivity and concentrations
    Output::ptr _solution;

    // private methods
    const double computeOmega(const precIndex j, const timeIndex n) const;

    const double computeZetaHat(const precIndex j, const timeIndex n) const;

    const double computePower(const timeIndex n, const double alpha) const;

    const double computeA1(const timeIndex n) const;

    const double computeB1(const timeIndex n) const;

    const double computeDT(const timeIndex n) const;

    const double computeGamma(const timeIndex n) const;

    const double getPower(const timeIndex n) const {
      return _solution->getPower(n);
    }

    const double getConcentration(const precIndex j, const timeIndex n) const {
      return _solution->getConcentration(j, n);
    }

    const double getRho(const timeIndex n) const {
      return _solution->getRho(n);
    }

  public:
    Solver(Params::ptr parameters, Output::ptr solution);

    Params::ptr getParameters() const { return _params; }

    Output::ptr getSolution() const { return _solution; }

    // Reset solution
    void reset(Output::ptr solution) {
      _solution = solution;
      baseReset(solution);
    }

    // Advance one time step
    void step(const timeIndex n) override;

    // Solve the epke problem
    void solve() override;
  }; // class Solver
} // namespace epke

#endif
