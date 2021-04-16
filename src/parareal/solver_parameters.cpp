#include "solver_parameters.hpp"
#include "utility/interpolate.hpp"

const EPKEParameters EPKEParameters::interpolate(const timeBins& fine_time) const {

  precBins<Precursor::ptr> fine_precursors;
  
  for (const auto p : _precursors) {
    Precursor::ptr fine_precursor
      = std::make_shared<Precursor>(util::interpolate(_time,
						      p->decayConstant(),
						      fine_time),
				    util::interpolate(_time,
						      p->delayedFraction(),
						      fine_time));
    fine_precursors.push_back(fine_precursor);
  }
  
  return EPKEParameters(fine_time,
			util::interpolate(_time, _rho_imp,  fine_time),
			util::interpolate(_time, _gen_time, fine_time),
			util::interpolate(_time, _pow_norm, fine_time),
			util::interpolate(_time, _beta_eff, fine_time),
			util::interpolate(_time, _lambda_h, fine_time),
			fine_precursors);
}
