#ifndef _EPKE_INPUT_HEADER_
#define _EPKE_INPUT_HEADER_

#include "precursor.hpp"
#include "pugixml.hpp"
#include "solver.hpp"

#include <string>

class Input {
private:
  std::string input_file_name;
  Solver::TimeIndex n_steps;

  std::vector<double> const loadVectorData(const pugi::xml_node& node);

  std::vector<Precursor::ptr> const loadPrecursors(
      const pugi::xml_node& precursors_node, const Solver::TimeIndex n_steps);

  const Solver::precBins<Solver::timeBins> loadConcentrationHistories(
      const pugi::xml_node& concentrations_node);

public:
  Input(std::string input_file_name) : input_file_name(input_file_name) {}

  void execute();
};

#endif