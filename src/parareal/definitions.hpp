#ifndef _PARAREAL_DEFINITIONS_HEADER_
#define _PARAREAL_DEFINITIONS_HEADER_

#include <vector>

namespace para {
  template <typename T>
  using precBins  = std::vector<T>;      // binning over precursor groups
  using timeBins  = std::vector<double>; // binning over time variable
  using timeIndex = uint32_t;
  using precIndex = uint8_t;
} // namespace para

#endif