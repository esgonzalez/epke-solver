#include "input.hpp"
#include "pugixml.hpp"
#include "solver.hpp"

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

void Input::execute() {
  pugi::xml_document input_file;
  pugi::xml_parse_result load_result
    = input_file.load_file(input_file_name.c_str());

  if (!load_result) {
    std::cout << load_result.description() << std::endl;
    throw;
  }

  std::cout << "Reading input file: " << input_file_name << std::endl;

  pugi::xml_node epke_node = input_file.child("epke_input");
  pugi::xml_node precursors_node = epke_node.child("precursors");
  pugi::xml_node time_node = epke_node.child("time");
  pugi::xml_node gen_time_node = epke_node.child("gen_time");
  pugi::xml_node pow_norm_node = epke_node.child("pow_norm");
  pugi::xml_node rho_imp_node = epke_node.child("rho_imp");
  pugi::xml_node beta_eff_node = epke_node.child("beta_eff");
  pugi::xml_node lambda_h_node = epke_node.child("lambda_h");

  timeBins<double> time = loadVectorData(time_node);
  timeBins<double> gen_time = loadVectorData(gen_time_node);
  timeBins<double> pow_norm = loadVectorData(pow_norm_node);
  timeBins<double> rho_imp = loadVectorData(rho_imp_node);
  timeBins<double> beta_eff = loadVectorData(beta_eff_node);
  timeBins<double> lambda_h = loadVectorData(lambda_h_node);
  precBins<Precursor::ptr> precursors
    = loadPrecursors(precursors_node, time.size());

  Solver solver(
    time, gen_time, pow_norm, rho_imp, beta_eff, lambda_h, precursors);

  double theta = epke_node.attribute("theta").as_double();
  double gamma_d = epke_node.attribute("gamma_d").as_double();
  double init_pow = epke_node.attribute("initial_power").as_double();
  double eta = epke_node.attribute("eta").as_double();

  solver.solve(theta, gamma_d, init_pow, eta);

  // build the xml document
  std::string outpath = "examples/epke_output.xml";
  std::ofstream out(outpath);
  pugi::xml_document doc;
  solver.buildXMLDoc( doc );
  std::cout << "writing output to " << outpath << std::endl;
  doc.save( out );
}

const std::vector<double> Input::loadVectorData( const pugi::xml_node &node ) {
  std::vector<double> result;
  std::istringstream iss(node.text().get());
  for (double s; iss >> s;) {
    result.push_back(s);
  }
  return result;
}

std::vector<Precursor::ptr> const Input::loadPrecursors(
    const pugi::xml_node &precursors_node, const uint32_t n_steps ) {
  std::vector<Precursor::ptr> precursors;
  for (auto precursor_node : precursors_node.children() ) {
    timeBins<double> decay_constant = std::vector<double>(
      n_steps, precursor_node.attribute("decay_constant").as_double());
    timeBins<double> beta = std::vector<double>(
      n_steps, precursor_node.attribute("beta").as_double());
    precursors.push_back(std::make_shared<Precursor>(decay_constant, beta));
  }
  return precursors;
}