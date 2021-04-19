#include "solver_parameters.hpp"

#include "utility/interpolate.hpp"
#include "utility/load_data.hpp"

EPKEParameters::EPKEParameters(const pugi::xml_node& epke_node) :
  _time(util::loadVectorData(epke_node.child("time"),
			     epke_node.attribute("n_steps").as_int())),
  _rho_imp(util::loadVectorData(epke_node.child("rho_imp"),
				epke_node.attribute("n_steps").as_int())),
  _gen_time(util::loadVectorData(epke_node.child("gen_time"),
				 epke_node.attribute("n_steps").as_int())),
  _pow_norm(util::loadVectorData(epke_node.child("pow_norm"),
				 epke_node.attribute("n_steps").as_int())),
  _beta_eff(util::loadVectorData(epke_node.child("beta_eff"),
				 epke_node.attribute("n_steps").as_int())),
  _lambda_h(util::loadVectorData(epke_node.child("lambda_h"),
				 epke_node.attribute("n_steps").as_int())),
  _p_history(util::loadVectorData(epke_node.child("p_history"),
				  epke_node.attribute("n_steps").as_int())),
  _zeta_histories(util::loadZetaHistories(epke_node.child("zeta_histories"))),
  _theta(epke_node.attribute("theta").as_double()),
  _gamma_d(epke_node.attribute("gamma_d").as_double()),
  _eta(epke_node.attribute("eta").as_double())
{
  timeIndex n_steps = epke_node.attribute("n_steps").as_int();

  // create the vector of precursors
  for (const auto& prec_node : epke_node.child("precursors")) {
    _precursors.push_back(std::make_shared<Precursor>(prec_node, n_steps));
  }
}

const EPKEParameters EPKEParameters::interpolate(const timeBins& fine_time) const {

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
			util::interpolate(_time, _rho_imp,  fine_time),
			util::interpolate(_time, _gen_time, fine_time),
			util::interpolate(_time, _pow_norm, fine_time),
			util::interpolate(_time, _beta_eff, fine_time),
			util::interpolate(_time, _lambda_h, fine_time),
			_p_history,
			_zeta_histories,
			fine_precursors,
			_theta,
			_gamma_d,
			_eta);
}
