#ifndef _PARAREAL_INPUT_HEADER_
#define _PARAREAL_INPUT_HEADER_

#include "parareal/precursor.hpp"
#include "pugi/pugixml.hpp"
#include "epke/epke_solver.hpp"

#include <string>

class Input {
private:
  std::string input_file_name;
  Solver::TimeIndex n_steps;

  std::vector<double> const loadVectorData(const pugi::xml_node& node);

  const Solver::precBins<Solver::timeBins> loadConcentrationHistories(
      const pugi::xml_node& concentrations_node);

public:
  Input(std::string input_file_name) : input_file_name(input_file_name) {}

  void execute();
};

#endif
