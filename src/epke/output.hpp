#ifndef _EPKE_OUTPUT_HEADER_
#define _EPKE_OUTPUT_HEADER_

#include <memory>

#include "parareal/solver_output.hpp"

namespace epke {

class EPKEOutput : public para::SolverOutput {
public:
  template <typename T>
  using precBins  = para::precBins<T>;
  using precIndex = para::precIndex;
  using ptr = std::shared_ptr<EPKEOutput>;

private:
  // time-dependent parameters
  timeBins           _power;          // reactor power
  timeBins           _pow_norm;       // Power normalization factor
  timeBins           _rho;            // reactivity with feedback
  precBins<timeBins> _concentrations; // precursor concentrations

public:
  EPKEOutput(const precIndex n_precursors,
	     const timeIndex n_time_steps,
	     const timeIndex n_start,
	     const timeIndex n_stop)
    : SolverOutput(n_time_steps, n_start, n_stop),
      _power(n_time_steps, 0.),
      _pow_norm(n_time_steps, 0.),
      _rho(n_time_steps, 0.),
      _concentrations(n_precursors, timeBins(n_time_steps, 0.)) {}

  // Construct from data vectors
  EPKEOutput(const timeIndex           n_start,
	     const timeIndex           n_stop,
	     const timeBins&           time,
	     const precBins<timeBins>& concentrations,
	     const timeBins&           power,
	     const timeBins&           pow_norm,
	     const timeBins&           rho)
    : SolverOutput(time, n_start, n_stop),
      _power(power),
      _pow_norm(pow_norm),
      _rho(rho),
      _concentrations(concentrations) {}

  void setPower(const timeIndex n, const double val) { _power[n] = val; }
  void setPowNorm(const timeIndex n, const double val) { _pow_norm[n] = val; }
  void setRho(const timeIndex n, const double val) { _rho[n] = val; }
  void setConcentration(const precIndex k, const timeIndex n, const double val)
  { _concentrations[k][n] = val; }

  const double getPower(const timeIndex n) const { return _power.at(n); }
  const double getPowNorm(const timeIndex n) const { return _pow_norm.at(n); }
  const double getRho(const timeIndex n)   const { return _rho.at(n);   }
  const double getConcentration(const precIndex k, const timeIndex n) const {
    return _concentrations.at(k).at(n);
  }

  const timeBins& getPower() const { return _power; }
  const timeBins& getPowNorm() const { return _pow_norm; }
  const timeBins& getRho() const { return _rho; }
  const precBins<timeBins>& getConcentrations() const {
    return _concentrations;
  }

  const precIndex getNumPrecursors() const { return _concentrations.size(); }

  SolverOutput::ptr createPrecomputedImpl(const timeIndex n,
					  const timeIndex n_fine_per_coarse)
    const override;

  EPKEOutput::ptr coarsen(const timeBins& coarse_time) const;

  void resize(const timeIndex n_steps) override;

  void updateCoarseImpl(const timeIndex n, SolverOutput::ptr coarse) override;

  void updatePararealImpl(const timeIndex n,
			  SolverOutput::ptr new_coarse,
			  SolverOutput::ptr fine_coarsened,
			  SolverOutput::ptr old_coarse) override;

  void writeToXML(pugi::xml_document& doc) const override;
};

} // namespace epke

#endif
