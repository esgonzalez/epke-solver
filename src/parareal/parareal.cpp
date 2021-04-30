#include <vector>
#include <cassert>

#include "parareal/parareal.hpp"

#include "utility/interpolate.hpp"

para::Parareal::Parareal(const pugi::xml_node& parareal_node) :
  _solver(parareal_node.child("epke_input"),
	  parareal_node.child("epke_output")),
  _outpath(parareal_node.attribute("outpath").value()),
  _max_iterations(parareal_node.attribute("max_iterations").as_int()),
  _n_fine_per_coarse(parareal_node.attribute("n_fine_per_coarse").as_int()) {}

para::timeBins para::Parareal::generateFineTime(para::timeIndex const n) {
  timeBins fine_time;

  auto params = _solver.getParameters();

  assert(n < params->getNumTimeSteps() - 1);

  if (params->getInterpolated()) {
    return params->getTime();
  }

  // push back the time steps before index n
  for (int i = 0; i < n; i++) {
    fine_time.push_back(params->getTime(i));
  }

  // push back the new time steps
  for (const auto& t : util::linspace(params->getTime(n),
				      params->getTime(n+1),
				      _n_fine_per_coarse+1)) {
    fine_time.push_back(t);
  }

  return fine_time;
}

void para::Parareal::solve() {
  auto params = _solver.getParameters();

  timeIndex coarse_size =
    params->getInterpolated() ? 1 : params->getNumTimeSteps() - 1;

  // loop over each index of the precomputed values
  for (timeIndex n = 0; n < coarse_size; n++) {
    const auto fine_time   = generateFineTime(n);
    const auto fine_solver = _solver.createFineSolver(fine_time, n);
    const auto fine_output = fine_solver->solve();
  }

  // have the coarse solver assemble the global output
  _global_output = _solver.assembleGlobalOutput();
}


void para::Parareal::writeToXML(pugi::xml_document& doc) const {
  std::ofstream out(_outpath);

  // create the root level node
  doc.append_child("parareal");

  auto params = _solver.getParameters();

  if (params->getInterpolated()) {
    _global_output->writeToXML(doc);
    params->writeToXML(doc);
  }
  else {
    auto global_coarsened = coarsen(_global_output, params->getTime());
    global_coarsened->writeToXML(doc);
 }

  doc.save(out);
}
