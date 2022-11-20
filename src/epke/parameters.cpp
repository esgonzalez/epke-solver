#include <iomanip>
#include <sstream>

#include "epke/parameters.hpp"

#include "utility/interpolate.hpp"
#include "pugi/pugixml.hpp"

para::SolverParameters::ptr
epke::EPKEParameters::interpolateImpl(const timeBins& fine_time) const {
  using util::interpolate;

  // interpolate the precursors
  precBins<Precursor::ptr> fine_precursors;

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
				     _eta);
}

void epke::EPKEParameters::writeToXML(pugi::xml_document& doc) const {
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
