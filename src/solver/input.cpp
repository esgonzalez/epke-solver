#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "input.hpp"

#include "parareal.hpp"
#include "solver_parameters.hpp"

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

  timeBins<double> time = loadVectorData(epke_node.child("time"));
  timeBins<double> rho_imp = loadVectorData(epke_node.child("rho_imp"));
  timeBins<double> gen_time = loadVectorData(epke_node.child("gen_time"));
  timeBins<double> pow_norm = loadVectorData(epke_node.child("pow_norm"));
  timeBins<double> beta_eff = loadVectorData(epke_node.child("beta_eff"));
  timeBins<double> lambda_h = loadVectorData(epke_node.child("lambda_h"));
  timeBins<double> p_history = loadVectorData(epke_node.child("p_history"));

  Solver::precBins<Precursor::ptr> precursors =
      loadPrecursors(epke_node.child("precursors"), time.size());

  Solver::precBins<Solver::timeBins> concentration_histories = loadConcentrationHistories(
      epke_node.child("concentration_histories"));

  EPKEParameters epke_parameters(time, gen_time, pow_norm, rho_imp,
				 beta_eff, lambda_h, precursors);
  
  Solver solver(time, gen_time, pow_norm, rho_imp,
		beta_eff, lambda_h, precursors);

  // Set the precomputed values if the simulation doesn't start at t=0
  solver.setPrecomputedValues(p_history, concentration_histories);

  // Grab additional simulation parameters
  double theta = epke_node.attribute("theta").as_double();
  double gamma_d = epke_node.attribute("gamma_d").as_double();
  double eta = epke_node.attribute("eta").as_double();

  // Run the EPKE solver
  solver.solve(theta, gamma_d, eta);

  // build the xml document
  std::string outpath = epke_node.attribute("outpath").value();
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

const Solver::precBins<Solver::timeBins> Input::loadConcentrationHistories(
    const pugi::xml_node& concentrations_node) {

  Solver::precBins<Solver::timeBins> concentration_histories;

  for (auto history_node : concentrations_node.children()) {
    std::vector<double> concentration_history;

    std::istringstream iss(history_node.text().get());
    for (double s; iss >> s;) {
     concentration_history.push_back(s);
    }

    concentration_histories.push_back(concentration_history);
  }

  return concentration_histories;
}
