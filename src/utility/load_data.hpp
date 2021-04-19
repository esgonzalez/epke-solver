#ifndef _UTILITY_LOAD_DATA_HEADER_
#define _UTILITY_LOAD_DATA_HEADER_

#include <sstream>

#include "pugi/pugixml.hpp"
#include "parareal/definitions.hpp"

namespace util {

inline const para::timeBins
loadVectorData(const pugi::xml_node& node, para::timeIndex n_steps) {
  para::timeBins result;

  // check to see if parameter is constant in time
  double value = node.attribute("value").as_double();

  if (value) {
    result = para::timeBins(n_steps, value);
  } else {
    std::istringstream iss(node.text().get());
    for (double s; iss >> s;) {
      result.push_back(s);
    }
  }

  return result;
}

} // namespace util

#endif
