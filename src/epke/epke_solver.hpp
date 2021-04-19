#ifndef _EPKE_EPKE_SOLVER_HEADER_
#define _EPKE_EPKE_SOLVER_HEADER_

#include "parareal/precursor.hpp"
#include "parareal/solver_parameters.hpp"
#include "parareal/definitions.hpp"
#include "pugi/pugixml.hpp"

class Solver {
public:
  template <typename T>
  using precBins  = para::precBins<T>;      // binning over precursor groups
  using timeBins  = para::timeBins; // binning over time variable
  using timeIndex = para::timeIndex;
  using precIndex = para::precIndex;

private:
  // input member variables
  const EPKEParameters params;

  // output member variables
  timeBins power; // power at various points in time
  timeBins rho;   // reactivity with feedback
  precBins<timeBins> concentrations; // precursor concentrations

  // TODO: move this to EPKEParameters
  timeBins delta_t;

  // time-dependent local variables (rewritten to every time step)
  precBins<double> omega, zeta_hat;

  // private methods
  const double computeOmega(const precIndex k,
			    const timeIndex n,
			    const double    w,
			    const double    gamma) const;

  const double computeZetaHat(const precIndex k,
			      const timeIndex n,
			      const double    w,
			      const double    gamma) const;

  const double computePower(const timeIndex n,
			    const double    alpha,
			    const double    gamma);

  const double computeABC(const timeIndex                  n,
			  const double                     alpha,
			  const std::pair<double, double>& a1b1,
			  const double                     tau,
			  const double                     s_hat_d,
			  const double                     s_d_prev) const;

  std::pair<double, double> computeA1B1(const timeIndex n, const double gamma);

  const bool acceptTransformation(const timeIndex n,
				  const double    alpha,
				  const double    gamma) const;

public:
  Solver(const EPKEParameters parameters);

  void solve();

  void buildXMLDoc(pugi::xml_document& doc) const;
};

#endif
