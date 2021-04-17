#ifndef _EPKE_EPKE_SOLVER_HEADER_
#define _EPKE_EPKE_SOLVER_HEADER_

#include "parareal/precursor.hpp"
#include "parareal/solver_parameters.hpp"
#include "pugi/pugixml.hpp"

#include <vector>

class Solver {
public:
  template <typename T>
  using precBins  = std::vector<T>;      // binning over precursor groups
  using timeBins  = std::vector<double>; // binning over time variable
  using TimeIndex = uint32_t;
  using PrecIndex = uint8_t;

private:
  // input member variables
  const EPKEParameters params;
  
  // output member variables
  timeBins power; // power at various points in time
  timeBins rho;   // reactivity with feedback (to be calculated)
  timeBins rho_d; // difference between imposed and calculated reactivities

  // TODO: move this to EPKEParameters
  timeBins delta_t;

  // power and concentration histories provided at start of simulation
  timeBins p_history;
  precBins<timeBins> concentrations; // precursor concentrations

  // time-dependent local variables (rewritten to every time step)
  precBins<double> omega, zeta_hat;

  // private methods
  const double computeOmega(
      const PrecIndex k, const TimeIndex n,
      const double w, const double gamma) const;
  const double computeZetaHat(
      const PrecIndex& k, const TimeIndex& n, const double& w,
      const double& gamma) const;
  const double computePower(
      const TimeIndex& n, const double& theta, const double& gamma_d,
      const double& eta, const double& alpha, const double& gamma);
  const double computeABC(
      const TimeIndex& n, const double theta, const double& alpha,
      const std::pair<double, double>& a1b1, const double& tau,
      const double& s_hat_d, const double& s_d_prev) const;
  std::pair<double, double> computeA1B1(
      const TimeIndex& n, const double gamma_d, const double eta,
      const double& gamma);
  const bool acceptTransformation(
      const TimeIndex& n, const double& alpha, const double& gamma) const;

public:
  Solver(const EPKEParameters parameters);

  // Set values that have been precomputed (for the parareal solver) so the
  // simulation does not have to begin at t=0
  void setPrecomputedValues(const timeBins& p_history,
			    const precBins<timeBins>& concentration_histories);

  // TODO: Move theta, gamma_d, and eta to EPKEParameters
  void solve(const double theta, const double gamma_d, const double eta);
  
  void buildXMLDoc(pugi::xml_document& doc) const;
};

#endif
