#include <algorithm>
#include <iterator>
#include <iomanip>
#include <sstream>

#include "epke/output.hpp"

#include "utility/interpolate.hpp"
#include "pugi/pugixml.hpp"

para::SolverOutput::ptr
epke::EPKEOutput::createPrecomputedImpl(const timeIndex n,
					const timeIndex n_fine_per_coarse) const
{
  timeIndex n_start = n * n_fine_per_coarse + 1;
  timeIndex n_stop = (n+1) * n_fine_per_coarse + 1;

  timeBins fine_time = util::linspace(_time.front(), _time.at(n), n_start);
  timeBins fine_power = util::interpolate(_time, _power, fine_time);
  timeBins fine_pow_norm = util::interpolate(_time, _pow_norm, fine_time);
  timeBins fine_rho   = util::interpolate(_time, _rho, fine_time);

  fine_time.resize(n_stop);
  fine_power.resize(n_stop);
  fine_pow_norm.resize(n_stop);
  fine_rho.resize(n_stop);

  precBins<timeBins> fine_concentrations;

  for (precIndex k = 0; k < getNumPrecursors(); k++) {
    timeBins fine_concentration = util::interpolate(_time,
						    _concentrations[k],
						    fine_time);
    fine_concentration.resize(n_stop);
    fine_concentrations.push_back(fine_concentration);
  }

  return std::make_shared<EPKEOutput>(n_start,
				      n_stop,
				      fine_time,
				      fine_concentrations,
				      fine_power,
				      fine_pow_norm,
				      fine_rho);
}

epke::EPKEOutput::ptr
epke::EPKEOutput::coarsen(const timeBins& coarse_time) const {
  timeBins coarse_power(coarse_time.size());
  timeBins coarse_pow_norm(coarse_time.size());
  timeBins coarse_rho(coarse_time.size());
  precBins<timeBins> coarse_concentrations(getNumPrecursors(),
					   timeBins(coarse_time.size()));

  for (timeIndex n_coarse = 0; n_coarse < coarse_time.size(); n_coarse++) {
    auto it = std::find(_time.cbegin(), _time.cend(), coarse_time.at(n_coarse));

    timeIndex n_fine = std::distance(_time.cbegin(), it);

    coarse_power[n_coarse] = _power.at(n_fine);
    coarse_pow_norm[n_coarse] = _pow_norm.at(n_fine);
    coarse_rho[n_coarse] = _rho.at(n_fine);

    for (precIndex k = 0; k < getNumPrecursors(); k++) {
      coarse_concentrations[k][n_coarse] = _concentrations.at(k).at(n_fine);
    }
  }

  return std::make_shared<EPKEOutput>(0,
				      coarse_time.size(),
				      coarse_time,
				      coarse_concentrations,
				      coarse_power,
				      coarse_pow_norm,
				      coarse_rho);
}

void epke::EPKEOutput::resize(const timeIndex n_steps) {
  _time.resize(n_steps);
  _power.resize(n_steps);
  _pow_norm.resize(n_steps);
  _rho.resize(n_steps);

  for (precIndex k = 0; k < getNumPrecursors(); k++) {
    _concentrations[k].resize(n_steps);
  }
}

void epke::EPKEOutput::updateCoarseImpl(const timeIndex n,
					para::SolverOutput::ptr coarse) {
  EPKEOutput::ptr coarse_output = std::static_pointer_cast<EPKEOutput>(coarse);

  coarse_output->setTime(n, getTime(n));
  coarse_output->setPower(n, getPower(n));
  coarse_output->setPowNorm(n, getPowNorm(n));
  coarse_output->setRho(n, getRho(n));
  for (precIndex j = 0; j < getNumPrecursors(); j++) {
    coarse_output->setConcentration(j, n, getConcentration(j,n));
  }
}

void epke::EPKEOutput::updatePararealImpl(const timeIndex n,
					  para::SolverOutput::ptr new_coarse,
					  para::SolverOutput::ptr fine_coarsened,
					  para::SolverOutput::ptr old_coarse) {

  auto new_co = std::static_pointer_cast<EPKEOutput>(new_coarse);
  auto fine_co = std::static_pointer_cast<EPKEOutput>(fine_coarsened);
  auto old_co = std::static_pointer_cast<EPKEOutput>(old_coarse);

  setPower(n,
	   new_co->getPower(n) +
	   fine_co->getPower(n) -
	   old_co->getPower(n));

  setPowNorm(n,
	     new_co->getPowNorm(n) +
	     fine_co->getPowNorm(n) -
	     old_co->getPowNorm(n));

  for (precIndex j = 0; j < getNumPrecursors(); j++) {
    setConcentration(j, n,
		     new_co->getConcentration(j,n) +
		     fine_co->getConcentration(j,n) -
		     old_co->getConcentration(j,n));
  }

  setRho(n,
	 new_co->getRho(n) +
	 fine_co->getRho(n) -
	 old_co->getRho(n));
}

void epke::EPKEOutput::writeToXML(pugi::xml_document& doc) const {
  pugi::xml_node parareal_node = doc.child("parareal");
  pugi::xml_node output_node   = parareal_node.append_child("epke_output");
  pugi::xml_node time_node     = output_node.append_child("time");
  pugi::xml_node power_node    = output_node.append_child("power");
  pugi::xml_node pow_norm_node = output_node.append_child("pow_norm");
  pugi::xml_node rho_node      = output_node.append_child("rho");
  pugi::xml_node concs_node    = output_node.append_child("concentrations");

  std::ostringstream time_str, power_str, pow_norm_str, rho_str, conc_str;

  parareal_node.append_attribute("solve_time") = _solve_time;

  for (int n = 0; n < getNumTimeSteps(); n++) {
    time_str     << std::setprecision(6)  << _time.at(n);
    power_str    << std::setprecision(12) << _power.at(n);
    pow_norm_str << std::setprecision(12) << _pow_norm.at(n);
    rho_str      << std::setprecision(12) << _rho.at(n);

    if (n != getNumTimeSteps() - 1) {
      time_str     << " ";
      power_str    << " ";
      pow_norm_str << " ";
      rho_str      << " ";
    }
  }

  time_node.text()     = time_str.str().c_str();
  power_node.text()    = power_str.str().c_str();
  pow_norm_node.text() = pow_norm_str.str().c_str();
  rho_node.text()      = rho_str.str().c_str();

  for (int k = 0; k < getNumPrecursors(); k++) {
    conc_str.str("");
    conc_str.clear();
    pugi::xml_node conc_node = concs_node.append_child("concentration");
    conc_node.append_attribute("k") = k;
    for (int n = 0; n < getNumTimeSteps(); n++) {
      conc_str << std::setprecision(12) << _concentrations.at(k).at(n);

      if (n != getNumTimeSteps() - 1) {
        conc_str << " ";
      }
    }
    conc_node.text() = conc_str.str().c_str();
  }
}
