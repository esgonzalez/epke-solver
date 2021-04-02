#ifndef _EPKE_SOLVER_HEADER_
#define _EPKE_SOLVER_HEADER_

#include "precursor.hpp"
#include "pugixml.hpp"

#include <vector>

class Solver {
public:
  template <typename T>
  using precBins = std::vector<T>;      // binning over precursor groups
  using timeBins = std::vector<double>; // binning over time variable
  using TimeIndex = uint32_t;
  using PrecIndex = uint8_t;

private:
  // input member variables
  const timeBins time;     // vector of time points
  const timeBins gen_time; // average neutron generation time (Lambda)
  const timeBins pow_norm; // power normalization factor (f_fp)
  const timeBins rho_imp;  // the imposed reactivity (without feedback)
  const timeBins beta_eff; // the total delayed neutron fraction
  const timeBins lambda_h;
  const precBins<Precursor::ptr> precursors;

  // output member variables
  timeBins power; // power at various points in time
  timeBins rho;   // reactivity with feedback (to be calculated)
  timeBins rho_d;
  timeBins delta_t;
  timeBins p_history;
  precBins<timeBins> concentrations; // precursor concentrations

  // time-dependent local variables (rewritten to every time step)
  precBins<double> omega, zeta_hat;

  // private methods
  const double computeOmega(
      const PrecIndex& k, const TimeIndex& n, const double& w,
      const double& gamma) const;
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
  Solver(
      const timeBins& time, const timeBins& gen_time, const timeBins& pow_norm,
      const timeBins& rho_imp, const timeBins& beta_eff,
      const timeBins& lambda_h, const timeBins& p_history,
      const precBins<Precursor::ptr>& precursors,
      const precBins<timeBins>& concentration_histories);

  void solve(const double theta, const double gamma_d, const double eta);

  void buildXMLDoc(pugi::xml_document& doc) const;
};

#endif