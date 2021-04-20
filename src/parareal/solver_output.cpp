#include "solver_output.hpp"

#include "utility/interpolate.hpp"
#include "utility/load_data.hpp"

epke::EPKEOutput::EPKEOutput(const pugi::xml_node& output_node) :
  _power(util::loadVectorData(output_node.child("power"))),
  _rho(util::loadVectorData(output_node.child("rho"))),
  _concentrations(util::loadZetas(output_node.child("concentrations"))) {}

const epke::EPKEOutput
epke::EPKEOutput::createPrecomputed(const timeIndex n) const {

  const timeBins power(_power.begin(), _power.begin() + n);
  const timeBins rho(_rho.begin(), _rho.begin() + n);
  precBins<timeBins> concentrations;
  for (precIndex k = 0; k < getNumPrecursors(); k++) {
    const timeBins concentration(_concentrations[k].begin(),
				 _concentrations[k].begin() + n);
    concentrations.push_back(concentration);
  }

  return EPKEOutput(power, rho, concentrations);
}
