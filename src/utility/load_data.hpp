#ifndef _UTILITY_LOAD_DATA_HEADER_
#define _UTILITY_LOAD_DATA_HEADER_

#include <sstream>

#include "pugi/pugixml.hpp"
#include "parareal/definitions.hpp"

namespace util {

inline const para::timeBins
loadVectorData(const pugi::xml_node& node, const para::timeIndex n_steps) {
  para::timeBins result;

  // check to see if parameter is constant in time
  double value = node.attribute("value").as_double();

  if (value) {
    result = para::timeBins(n_steps, value);
  } else {
    result = para::timeBins(n_steps, 0.);
    std::istringstream iss(node.text().get());
    para::timeIndex n = 0;
    for (double s; iss >> s;) {
      result[n] = s;
      n++;
    }
  }

  return result;
}

inline const para::timeBins
loadVectorData(const pugi::xml_node& node) {
  para::timeBins result;

  std::istringstream iss(node.text().get());
  for (double s; iss >> s;) {
    result.push_back(s);
  }

  return result;
}

inline const para::precBins<para::timeBins>
loadZetas(const pugi::xml_node& histories_node, const para::timeIndex n_steps) {

  para::precBins<para::timeBins> concentration_histories;

  for (auto history_node : histories_node.children()) {
    para::timeBins concentration_history(n_steps, 0.);

    para::timeIndex n = 0;
    std::istringstream iss(history_node.text().get());
    for (double s; iss >> s;) {
      concentration_history[n] = s;
      n++;
    }

    concentration_histories.push_back(concentration_history);
  }

  return concentration_histories;
}

} // namespace util

#endif
