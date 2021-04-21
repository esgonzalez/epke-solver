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

inline const para::timeBins
loadVectorData(const pugi::xml_node& node) {
  para::timeBins result;

  std::istringstream iss(node.text().get());
  for (double s; iss >> s;) {
    result.push_back(s);
  }

  return result;
}

inline const para::precBins<para::timeBins> loadZetas(
    const pugi::xml_node& histories_node) {

  template <typename T> using timeBins = std::vector<T>;
  template <typename T> using precBins = std::vector<T>;

  timeBins<precBins<double>> concentration_histories;

  for (auto history_node : histories_node.children()) {
    std::vector<double> concentration_history;

    std::istringstream iss(history_node.text().get());
    for (double s; iss >> s;) {
     concentration_history.push_back(s);
    }

    concentration_histories.push_back(concentration_history);
  }

  return concentration_histories;
}

} // namespace util

#endif
