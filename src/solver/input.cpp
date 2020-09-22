#include "input.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

void Input::execute() {
  pugi::xml_document input_file;
  pugi::xml_parse_result load_result =
      input_file.load_file(input_file_name.c_str());

  if (!load_result) {
    std::cout << load_result.description() << std::endl;
    throw;
  }

  std::cout << "Reading input file: " << input_file_name << std::endl;

  pugi::xml_node epke_node = input_file.child("epke_input");

  n_steps = epke_node.attribute("n_steps").as_int();

  timeBins<double> time     = loadVectorData(epke_node.child("time"));
  timeBins<double> rho_imp  = loadVectorData(epke_node.child("rho_imp"));
  timeBins<double> gen_time = loadVectorData(epke_node.child("gen_time"));
  timeBins<double> pow_norm = loadVectorData(epke_node.child("pow_norm"));
  timeBins<double> beta_eff = loadVectorData(epke_node.child("beta_eff"));
  timeBins<double> lambda_h = loadVectorData(epke_node.child("lambda_h"));

  Solver::precBins<Precursor::ptr> precursors =
      loadPrecursors(epke_node.child("precursors"), time.size());

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
  solver.buildXMLDoc(doc);
  std::cout << "writing output to " << outpath << std::endl;
  doc.save(out);
}

const std::vector<double> Input::loadVectorData(const pugi::xml_node& node) {
  std::vector<double> result;

  // check to see if parameter is constant in time
  double value = node.attribute("value").as_double();

  if (value) {
    result = std::vector<double>(n_steps, value);
  } else {
    std::istringstream iss(node.text().get());
    for (double s; iss >> s;) {
      result.push_back(s);
    }
  }

  return result;
}

std::vector<Precursor::ptr> const Input::loadPrecursors(
    const pugi::xml_node& precursors_node, const uint32_t n_steps) {
  std::vector<Precursor::ptr> precursors;
  for (auto precursor_node : precursors_node.children()) {
    timeBins<double> decay_constant = std::vector<double>(
        n_steps, precursor_node.attribute("decay_constant").as_double());
    timeBins<double> beta = std::vector<double>(
        n_steps, precursor_node.attribute("beta").as_double());
    precursors.push_back(std::make_shared<Precursor>(decay_constant, beta));
  }
  return precursors;
}