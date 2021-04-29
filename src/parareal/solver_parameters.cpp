#include <iomanip>

#include "solver_parameters.hpp"

#include "utility/interpolate.hpp"
#include "utility/load_data.hpp"

para::SolverParameters::SolverParameters(const pugi::xml_node& solver_node) :
  _time(util::loadVectorData(solver_node.child("time"),
			     solver_node.attribute("n_steps").as_int())),
  _precursors(buildPrecursors(solver_node.child("precursors"),
			      solver_node.attribute("n_steps").as_int())),
  _interpolated(solver_node.attribute("interpolated").as_bool()) {}

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

para::SolverParameters::ptr
EPKEParameters::interpolateImpl(const timeBins& fine_time) const {
  using util::interpolate;

  // interpolate the precursors
  precBins<Precursor::ptr> fine_precursors;

  if (_interpolated) {
    return std::make_shared<EPKEParameters>(*this);
  }

  for (const auto p : _precursors) {
    Precursor::ptr fine_precursor
      = std::make_shared<Precursor>(interpolate(_time,
						p->decayConstant(),
						fine_time),
				    interpolate(_time,
						p->delayedFraction(),
						fine_time));
    fine_precursors.push_back(fine_precursor);
  }

  return
    std::make_shared<EPKEParameters>(fine_time,
				     fine_precursors,
				     interpolate(_time, _rho_imp,  fine_time),
				     interpolate(_time, _gen_time, fine_time),
				     interpolate(_time, _pow_norm, fine_time),
				     interpolate(_time, _beta_eff, fine_time),
				     interpolate(_time, _lambda_h, fine_time),
				     _theta,
				     _gamma_d,
				     _eta,
				     true);
}

void EPKEParameters::writeToXML(pugi::xml_document& doc) const {
  pugi::xml_node parareal_node = doc.child("parareal");
  pugi::xml_node params_node   = parareal_node.append_child("epke_input");
  pugi::xml_node time_node     = params_node.append_child("time");
  pugi::xml_node rho_imp_node  = params_node.append_child("rho_imp");
  pugi::xml_node gen_time_node = params_node.append_child("gen_time");
  pugi::xml_node pow_norm_node = params_node.append_child("pow_norm");
  pugi::xml_node beta_eff_node = params_node.append_child("beta_eff");
  pugi::xml_node lambda_h_node = params_node.append_child("lambda_h");

  std::ostringstream time_str, rho_imp_str, gen_time_str, pow_norm_str,
    beta_eff_str, lambda_h_str;

  for (int n = 0; n < getNumTimeSteps(); n++) {
    time_str     << std::setprecision(6)  << _time.at(n);
    rho_imp_str  << std::setprecision(12) << _rho_imp.at(n);
    gen_time_str << std::setprecision(12) << _gen_time.at(n);
    pow_norm_str << std::setprecision(12) << _pow_norm.at(n);
    beta_eff_str << std::setprecision(12) << _beta_eff.at(n);
    lambda_h_str << std::setprecision(12) << _lambda_h.at(n);

    if (n != getNumTimeSteps() - 1) {
      time_str     << " ";
      rho_imp_str  << " ";
      gen_time_str << " ";
      pow_norm_str << " ";
      beta_eff_str << " ";
      lambda_h_str << " ";
    }
  }

  params_node.append_attribute("theta")   = _theta;
  params_node.append_attribute("gamma_d") = _gamma_d;
  params_node.append_attribute("eta")     = _eta;

  time_node.text()     = time_str.str().c_str();
  rho_imp_node.text()  = rho_imp_str.str().c_str();
  gen_time_node.text() = gen_time_str.str().c_str();
  pow_norm_node.text() = pow_norm_str.str().c_str();
  beta_eff_node.text() = beta_eff_str.str().c_str();
  lambda_h_node.text() = lambda_h_str.str().c_str();
}
