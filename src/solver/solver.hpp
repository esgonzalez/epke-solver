#ifndef _EPKE_SOLVER_HEADER_
#define _EPKE_SOLVER_HEADER_

#include "precursor.hpp"
#include "pugixml.hpp"

#include <vector>

template <typename T> using timeBins = std::vector<T>;
template <typename T> using precBins = std::vector<T>;

class Solver {
private:
  // input member variables
  const timeBins<double> time; // vector of time points we will be
  const timeBins<double> gen_time; // average neutron generation time (Lambda)
  const timeBins<double> pow_norm; // power normalization factor (f_fp)
  const timeBins<double> rho_imp; // the imposed reactivity (without feedback)
  const timeBins<double> beta_eff; // the total delayed neutron fractions
  const timeBins<double> lambda_h;
  const precBins<Precursor::ptr> precursors;
  precBins< timeBins<double> > concentrations; // precursor concentrations

  // output member variables
  timeBins<double> power; // power at various points in time
  timeBins<double> H; // I'm not sure what this is
  timeBins<double> rho; // with feedback (to be calculated)
  timeBins<double> rho_d;
  timeBins<double> delta_t;

  // these are rewritten to every time step
  precBins<double> omega, zeta_hat;
  double a1, b1, tau, s_hat_d, s_d_prev, alpha, w, gamma;

  double computeOmega( const uint8_t k, const uint32_t n ) const;
  double computeZetaHat( const uint8_t k, const uint32_t n ) const;
  double computePower(
    const uint32_t n, const double theta, const double gamma_d,
    const double eta );
  double computeABC( const uint32_t n, const double theta ) const;
  void computeA1B1( const uint32_t n, const double gamma_d, const double eta );
  bool acceptTransformation( const uint32_t n ) const;

public:
  Solver(
    const timeBins<double> time, const timeBins<double> gen_time,
    const timeBins<double> pow_norm, const timeBins<double> rho_imp,
    const timeBins<double> beta_eff, const timeBins<double> lambda_h,
    const precBins<Precursor::ptr> precursors );

  void solve(
    const double theta, const double gamma_d, const double init_pow,
    const double eta );
  void buildXMLDoc( pugi::xml_document& doc ) const;
};

#endif