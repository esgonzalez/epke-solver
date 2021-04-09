#include <vector>
#include <cassert>

#include "parareal/parareal.hpp"

#include "utility/interpolate.hpp"

Parareal::TimeBins Parareal::generateFineTime(Parareal::TimeIndex n) {
  TimeBins fine_time;

  assert(n < _coarse_time.size() - 1);

  // push back the time steps before index n
  for (int i = 0; i < n; i++) {
    fine_time.push_back(_coarse_time.at(i));
  }

  // push back the new time steps
  for (const auto& t : linspace(_coarse_time.at(n),
				_coarse_time.at(n+1),
				_n_fine_per_coarse+1)) {
    fine_time.push_back(t);
  }
    
  return fine_time;
}
