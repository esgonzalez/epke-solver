#include <iostream>
#include <vector>

#include "input.hpp"

#include "parareal/parareal.hpp"
#include "pugi/pugixml.hpp"

void Input::execute() {
  pugi::xml_document input_file;
  pugi::xml_parse_result load_result =
      input_file.load_file(input_file_name.c_str());

  if (!load_result) {
    std::cout << load_result.description() << std::endl;
    throw;
  }

  std::cout << "Reading input file: " << input_file_name << std::endl;

  para::Parareal parareal(input_file.child("parareal"));

  // Run the EPKE solver
  std::cout << "Solving..." << std::endl;
  parareal.solve();
  std::cout << "Completed solve." << std::endl;

  // build the xml document
  std::cout << "Writing output to " << parareal.getOutpath() << std::endl;
  parareal.writeToXML();
}
