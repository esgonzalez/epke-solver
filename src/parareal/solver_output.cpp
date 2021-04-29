#include <algorithm>
#include <iterator>
#include <iomanip>

#include "solver_output.hpp"

#include "utility/interpolate.hpp"
#include "utility/load_data.hpp"

para::SolverOutput::SolverOutput(const pugi::xml_node& precomp_node)
  : _time(util::loadVectorData(precomp_node.child("time"),
			       precomp_node.attribute("n_steps").as_int())),
    _concentrations(util::loadZetas(precomp_node.child("concentrations"))) {}

epke::EPKEOutput::EPKEOutput(const pugi::xml_node& precomp_node) :
  SolverOutput(precomp_node),
  _power(util::loadVectorData(precomp_node.child("power"))),
  _rho(util::loadVectorData(precomp_node.child("rho"))) {}

para::SolverOutput::ptr
epke::EPKEOutput::createPrecomputedImpl(const timeIndex n) const {
  const timeBins time(_time.begin(), _time.begin() + n + 1);
  const timeBins power(_power.begin(), _power.begin() + n + 1);
  const timeBins rho(_rho.begin(), _rho.begin() + n + 1);
  precBins<timeBins> concentrations;
  for (precIndex k = 0; k < getNumPrecursors(); k++) {
    const timeBins concentration(_concentrations[k].begin(),
				 _concentrations[k].begin() + n + 1);
    concentrations.push_back(concentration);
  }

  return std::make_shared<EPKEOutput>(time, concentrations, power, rho);
}

para::SolverOutput::ptr
epke::EPKEOutput::coarsenImpl(const timeBins& coarse_time) const {
  timeBins coarse_power(coarse_time.size());
  timeBins coarse_rho(coarse_time.size());
  precBins<timeBins> coarse_concentrations(getNumPrecursors(),
					   timeBins(coarse_time.size()));

  for (timeIndex n_coarse = 0; n_coarse < coarse_time.size(); n_coarse++) {
    auto it = std::find(_time.cbegin(), _time.cend(), coarse_time.at(n_coarse));

    timeIndex n_fine = std::distance(_time.cbegin(), it);

    coarse_power[n_coarse] = _power.at(n_fine);
    coarse_rho[n_coarse] = _rho.at(n_fine);

    for (precIndex k = 0; k < getNumPrecursors(); k++) {
      coarse_concentrations[k][n_coarse] = _concentrations.at(k).at(n_fine);
    }
  }

  return std::make_shared<EPKEOutput>(coarse_time,
				      coarse_concentrations,
				      coarse_power,
				      coarse_rho);
}

void epke::EPKEOutput::writeToXML(pugi::xml_document& doc) const {
  pugi::xml_node output_node = doc.append_child("epke_output");
  pugi::xml_node time_node = output_node.append_child("time");
  pugi::xml_node power_node = output_node.append_child("power");
  pugi::xml_node rho_node = output_node.append_child("rho");
  pugi::xml_node concs_node = output_node.append_child("concentrations");
  std::ostringstream time_str, power_str, rho_str, conc_str;

  for (int n = 0; n < getNumTimeSteps(); n++) {
    time_str << std::setprecision(6) << _time.at(n);
    power_str << std::setprecision(12) << _power.at(n);
    rho_str << std::setprecision(12) << _rho.at(n);
    if (n != getNumTimeSteps() - 1) {
      time_str << " ";
      power_str << " ";
      rho_str << " ";
    }
  }

  time_node.text() = time_str.str().c_str();
  power_node.text() = power_str.str().c_str();
  rho_node.text() = rho_str.str().c_str();

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
