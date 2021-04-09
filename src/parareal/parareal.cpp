#include <vector>

#include "parareal/parareal.hpp"

#include "utility/interpolate.hpp"

Parareal::TimeBins Parareal::generateFineTime(Parareal::TimeIndex n) {
  TimeBins fine_time;

  assert(n < _coarse_time.size() - 1);
  
  for (int i = 0; i < n; i++) {
    fine_time.push_back(_coarse_time.at(i));
  }

  TimeBins new_time = linspace(_coarse_time.at(n),
			       _coarse_time.at(n+1),
			       _n_fine_per_coarse+1);

  for (const auto& t : new_time) {
    fine_time.push_back(t);
  }
    
  return fine_time;
}
