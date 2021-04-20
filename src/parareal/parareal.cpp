#include <vector>
#include <cassert>

#include "parareal/parareal.hpp"

#include "utility/interpolate.hpp"

para::timeBins para::Parareal::generateFineTime(para::timeIndex const n) {
  timeBins fine_time;

  assert(n < _coarse_parameters.getNumTimeSteps() - 1);

  // push back the time steps before index n
  for (int i = 0; i < n; i++) {
    fine_time.push_back(_coarse_parameters.getTime(i));
  }

  // push back the new time steps
  for (const auto& t : util::linspace(_coarse_parameters.getTime(n),
				      _coarse_parameters.getTime(n+1),
				      _n_fine_per_coarse+1)) {
    fine_time.push_back(t);
  }

  return fine_time;
}
