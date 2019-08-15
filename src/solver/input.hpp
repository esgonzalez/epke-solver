#ifndef _EPKE_INPUT_HEADER_
#define _EPKE_INPUT_HEADER_

#include "precursor.hpp"
#include "pugixml.hpp"

#include <string>

class Input {
private:
  std::string input_file_name;

  std::vector<double> const loadVectorData( const pugi::xml_node &node );
  std::vector<Precursor::ptr> const loadPrecursors(
    const pugi::xml_node &precursors_node, const uint32_t n_steps );
public:
  Input(std::string input_file_name) : input_file_name(input_file_name) {}

  void execute();
};

#endif