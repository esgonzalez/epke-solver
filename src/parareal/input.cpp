#include <iostream>
#include <vector>

#include "input.hpp"

#include "parareal/parareal.hpp"
#include "epke/precursor.hpp"
#include "epke/parameters.hpp"
#include "epke/output.hpp"
#include "epke/solver.hpp"
#include "pugi/pugixml.hpp"
#include "utility/interpolate.hpp"
#include "utility/load_data.hpp"

void Input::execute() {
  pugi::xml_document input_file;
  pugi::xml_parse_result load_result =
      input_file.load_file(input_file_name.c_str());

  if (!load_result) {
    std::cout << load_result.description() << std::endl;
    throw;
  }

  std::cout << "Reading input file: " << input_file_name << std::endl;

  using Coarse = epke::Solver;
  using Fine   = epke::Solver;

  pugi::xml_node parareal_node = input_file.child("parareal");

  using util::loadVectorData;

  Coarse::Params::ptr coarse_params;
  Coarse::Output::ptr coarse_precomp;

  Fine::Params::ptr fine_params;
  Fine::Output::ptr fine_precomp;

  {
    // Create the coarse parameters from xml
    const pugi::xml_node params_node = parareal_node.child("epke_input");
    timeIndex n_steps = params_node.attribute("n_steps").as_int();

    timeBins time     = loadVectorData(params_node.child("time"), n_steps);
    timeBins rho_imp  = loadVectorData(params_node.child("rho_imp"), n_steps);
    timeBins gen_time = loadVectorData(params_node.child("gen_time"), n_steps);
    timeBins pow_norm = loadVectorData(params_node.child("pow_norm"), n_steps);
    timeBins beta_eff = loadVectorData(params_node.child("beta_eff"), n_steps);
    timeBins lambda_h = loadVectorData(params_node.child("lambda_h"), n_steps);
    double   theta    = params_node.attribute("theta").as_double();
    double   gamma_d  = params_node.attribute("gamma_d").as_double();
    double   eta      = params_node.attribute("eta").as_double();

    // Create the precursors
    const pugi::xml_node precs_node = params_node.child("precursors");
    precBins<epke::Precursor::ptr> precursors;
    timeBins lambda; timeBins beta;

    for (const auto& prec_node : precs_node) {
      lambda = loadVectorData(prec_node.child("decay_constant"), n_steps);
      beta   = loadVectorData(prec_node.child("delayed_fraction"), n_steps);
      precursors.push_back(std::make_shared<epke::Precursor>(lambda, beta));
    }

    // Create the coarse parameters object
    coarse_params = std::make_shared<Coarse::Params>(time,
						     precursors,
						     rho_imp,
						     gen_time,
						     pow_norm,
						     beta_eff,
						     lambda_h,
						     theta,
						     gamma_d,
						     eta);
  }

  {
    // Create the coarse solver initial conditions from xml
    const pugi::xml_node precomp_node = parareal_node.child("epke_output");
    timeIndex n_steps = precomp_node.attribute("n_steps").as_int();

    timeIndex n_start = precomp_node.attribute("n_start").as_int();
    timeIndex n_stop  = precomp_node.attribute("n_stop").as_int();
    timeBins time     = loadVectorData(precomp_node.child("time"), n_steps);
    timeBins power    = loadVectorData(precomp_node.child("power"), n_steps);
    timeBins pow_norm = loadVectorData(precomp_node.child("pow_norm"), n_steps);
    timeBins rho      = loadVectorData(precomp_node.child("rho"), n_steps);
    precBins<timeBins> concentrations =
      util::loadZetas(precomp_node.child("concentrations"), n_steps);

    coarse_precomp = std::make_shared<Coarse::Output>(n_start,
						      n_stop,
						      time,
						      concentrations,
						      power,
						      pow_norm,
						      rho);
  }

  Coarse::ptr coarse_solver =
    std::make_shared<Coarse>(coarse_params, coarse_precomp);

  // Create fine time
  timeIndex n_fine_per_coarse =
    parareal_node.attribute("n_fine_per_coarse").as_int();
  timeIndex n_fine =
    (coarse_solver->getNumTimeSteps() - 1) * n_fine_per_coarse + 1;
  timeBins fine_time =
    util::linspace(0., coarse_solver->getTime().back(), n_fine);

  {
    // Create the fine parameters
    fine_params = para::interpolate(coarse_params, fine_time);
  }

  {
    // Create the fine solver initial conditions from xml
    const pugi::xml_node precomp_node = parareal_node.child("epke_output");
    timeIndex n_steps = precomp_node.attribute("n_steps").as_int();

    timeIndex n_start = precomp_node.attribute("n_start").as_int();
    timeIndex n_stop  = precomp_node.attribute("n_stop").as_int();
    timeBins time     = loadVectorData(precomp_node.child("time"), n_steps);
    timeBins power    = loadVectorData(precomp_node.child("power"), n_steps);
    timeBins pow_norm = loadVectorData(precomp_node.child("pow_norm"), n_steps);
    timeBins rho      = loadVectorData(precomp_node.child("rho"), n_steps);
    precBins<timeBins> concentrations =
      util::loadZetas(precomp_node.child("concentrations"), n_steps);

    fine_precomp = std::make_shared<Fine::Output>(n_start,
						  n_stop,
						  time,
						  concentrations,
						  power,
						  pow_norm,
						  rho);

    fine_precomp->resize(n_fine);
  }

  Fine::ptr fine_solver = std::make_shared<Fine>(fine_params, fine_precomp);

  // Create the parareal solver
  para::Parareal<Coarse, Fine>
    parareal(coarse_solver,
	     fine_solver,
	     fine_precomp,
	     parareal_node.attribute("n_fine_per_coarse").as_int(),
	     parareal_node.attribute("max_iterations").as_int(),
	     parareal_node.attribute("outpath").value());

  // Run the EPKE solver
  std::cout << "Solving..." << std::endl;

  parareal.solve();
  std::cout << "Completed solve." << std::endl;

  // build the xml document
  std::cout << "Writing output to " << parareal.getOutpath() << std::endl;
  pugi::xml_document doc;
  parareal.writeToXML(doc);
}
