#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "input.hpp"

#include "epke/epke_solver.hpp"
#include "parareal/parareal.hpp"
#include "parareal/solver_parameters.hpp"
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

  pugi::xml_node parareal_node = input_file.child("parareal");
  pugi::xml_node epke_input_node = parareal_node.child("epke_input");
  pugi::xml_node epke_output_node = parareal_node.child("epke_output");
  EPKEParameters epke_params(epke_input_node);
  epke::EPKEOutput epke_precomputed(epke_output_node);
  epke::Solver epke_solver(epke_params, epke_precomputed);
  para::Parareal parareal(epke_params, epke_precomputed, epke_solver, 1, 1);

  // Run the EPKE solver
  std::cout << "Solving..." << std::endl;
  parareal.solve();
  std::cout << "Completed solve." << std::endl;

  // build the xml document
  std::string outpath = epke_input_node.attribute("outpath").value();
  std::ofstream out(outpath);
  pugi::xml_document doc;
  epke_solver.buildXMLDoc(doc);
  std::cout << "Writing output to " << outpath << std::endl;
  doc.save(out);
}
