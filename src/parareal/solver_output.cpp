#include <iomanip>

#include "solver_output.hpp"

#include "utility/interpolate.hpp"
#include "utility/load_data.hpp"

epke::EPKEOutput::EPKEOutput(const pugi::xml_node& output_node) :
  _power(util::loadVectorData(output_node.child("power"))),
  _rho(util::loadVectorData(output_node.child("rho"))),
  _concentrations(util::loadZetas(output_node.child("concentrations"))) {}

const epke::EPKEOutput
epke::EPKEOutput::createPrecomputed(const timeIndex n) const {

  const timeBins power(_power.begin(), _power.begin() + n + 1);
  const timeBins rho(_rho.begin(), _rho.begin() + n + 1);
  precBins<timeBins> concentrations;
  for (precIndex k = 0; k < getNumPrecursors(); k++) {
    const timeBins concentration(_concentrations[k].begin(),
				 _concentrations[k].begin() + n + 1);
    concentrations.push_back(concentration);
  }

  return EPKEOutput(power, rho, concentrations);
}

void epke::EPKEOutput::writeToXML(pugi::xml_document& doc) const {
  pugi::xml_node output_node = doc.append_child("epke_output");
  //pugi::xml_node time_node = output_node.append_child("time");
  pugi::xml_node power_node = output_node.append_child("power");
  pugi::xml_node rho_node = output_node.append_child("rho");
  pugi::xml_node concs_node = output_node.append_child("concentrations");
  //std::ostringstream time_str;
  std::ostringstream power_str, rho_str, conc_str;

  for (int n = 0; n < getNumTimeSteps(); n++) {
    //time_str << std::setprecision(6) << _time.at(n);
    power_str << std::setprecision(12) << _power.at(n);
    rho_str << std::setprecision(12) << _rho.at(n);
    if (n != getNumTimeSteps() - 1) {
      //time_str << " ";
      power_str << " ";
      rho_str << " ";
    }
  }

  //time_node.text() = time_str.str().c_str();
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
