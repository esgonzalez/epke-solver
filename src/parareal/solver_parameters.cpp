#include "solver_parameters.hpp"

#include "utility/interpolate.hpp"
#include "utility/load_data.hpp"

para::SolverParameters::SolverParameters(const pugi::xml_node& solver_node) :
  _time(util::loadVectorData(solver_node.child("time"),
			     solver_node.attribute("n_steps").as_int())),
  _precursors(buildPrecursors(solver_node.child("precursors"),
			      solver_node.attribute("n_steps").as_int())) {}

EPKEParameters::EPKEParameters(const pugi::xml_node& input_node) :
  SolverParameters(input_node),
  _rho_imp(util::loadVectorData(input_node.child("rho_imp"),
				input_node.attribute("n_steps").as_int())),
  _gen_time(util::loadVectorData(input_node.child("gen_time"),
				 input_node.attribute("n_steps").as_int())),
  _pow_norm(util::loadVectorData(input_node.child("pow_norm"),
				 input_node.attribute("n_steps").as_int())),
  _beta_eff(util::loadVectorData(input_node.child("beta_eff"),
				 input_node.attribute("n_steps").as_int())),
  _lambda_h(util::loadVectorData(input_node.child("lambda_h"),
				 input_node.attribute("n_steps").as_int())),
  _theta(input_node.attribute("theta").as_double()),
  _gamma_d(input_node.attribute("gamma_d").as_double()),
  _eta(input_node.attribute("eta").as_double()) {}

const EPKEParameters
EPKEParameters::interpolate(const timeBins& fine_time) const {
  // interpolate the precursors
  precBins<Precursor::ptr> fine_precursors;

  for (const auto p : _precursors) {
    Precursor::ptr fine_precursor
      = std::make_shared<Precursor>(util::interpolate(_time,
						      p->decayConstant(),
						      fine_time),
				    util::interpolate(_time,
						      p->delayedFraction(),
						      fine_time));
    fine_precursors.push_back(fine_precursor);
  }

  return EPKEParameters(fine_time,
			fine_precursors,
			util::interpolate(_time, _rho_imp,  fine_time),
			util::interpolate(_time, _gen_time, fine_time),
			util::interpolate(_time, _pow_norm, fine_time),
			util::interpolate(_time, _beta_eff, fine_time),
			util::interpolate(_time, _lambda_h, fine_time),
			_theta,
			_gamma_d,
			_eta);
}
