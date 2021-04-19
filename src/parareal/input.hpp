#ifndef _PARAREAL_INPUT_HEADER_
#define _PARAREAL_INPUT_HEADER_

#include "parareal/definitions.hpp"

#include <string>

class Input {
private:
  std::string input_file_name;

public:
  Input(std::string input_file_name) : input_file_name(input_file_name) {}

  void execute();
};

#endif
