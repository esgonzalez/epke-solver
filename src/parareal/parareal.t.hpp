#include <iomanip>
#include <sstream>
#include <fstream>
#include <chrono>

#include "utility/interpolate.hpp"
#include "pugi/pugixml.hpp"

using namespace para;

template <typename Coarse, typename Fine>
timeBins Parareal<Coarse, Fine>::generateFineTime(const timeIndex n) const {
  return util::linspace(_coarse_solver->getTime(0),
			_coarse_solver->getTime(n+1),
			(n+1) * _n_fine_per_coarse+1);
}

template <typename Coarse, typename Fine>
void Parareal<Coarse, Fine>::runCoarseSolver() {
  // Run the coarse solver
  _coarse_solver->solve();

  // Update the old coarse solution
  for (timeIndex n = 0; n < _coarse_solver->getNumTimeSteps(); n++) {
    updateCoarse(n, _coarse_solver->getSolution(), _old_coarse);
  }
}

template <typename Coarse, typename Fine>
void Parareal<Coarse, Fine>::runFineSolver(const timeIndex n) {
  timeBins fine_time = generateFineTime(n);

  // TODO: Move these to input.cpp and use global output to generate new precomp
  // Interpolate the precomputed values for the fine solver
  typename Output::ptr fine_precomp =
    createPrecomputed(_coarse_solver->getSolution(), n,	_n_fine_per_coarse);

  _fine_solver->reset(fine_precomp);
  _fine_solver->solve();

  assembleGlobalOutput();
}

template <typename Coarse, typename Fine>
void Parareal<Coarse, Fine>::update(const paraIndex k) {
  _fine_coarsened = _global_output->coarsen(_coarse_solver->getTime());

  double err = 0.;

  // set initial conditions for the power and reactivity vectors
  for (timeIndex n = 1; n < _coarse_solver->getNumTimeSteps(); n++) {
    _coarse_solver->step(n);

    // update new coarse solution
    updateCoarse(n, _coarse_solver->getSolution(), _new_coarse);

    // make parareal adjustments
    updateParareal(n,
		   _coarse_solver->getSolution(),
		   _new_coarse,
		   _fine_coarsened,
		   _old_coarse);

    // compute error for convergence
    //err += fabs(_new_coarse->getPower(n) -
    //_old_coarse->getPower(n)) / _coarse_solver.getNumTimeSteps();
  }

  //std::cout << "err: " << err << std::endl;

  // Swap the old and new coarse solutions
  std::swap(_new_coarse, _old_coarse);
}

template <typename Coarse, typename Fine>
void Parareal<Coarse, Fine>::assembleGlobalOutput() {
  for (timeIndex n_fine = _fine_solver->getStartTimeIndex();
       n_fine < _fine_solver->getStopTimeIndex(); n_fine++) {
    updateCoarse(n_fine, _fine_solver->getSolution(), _global_output);
  }
}

template <typename Coarse, typename Fine>
void Parareal<Coarse, Fine>::solve() {
  using timer = std::chrono::system_clock;
  using namespace std::chrono_literals;

  timer::time_point clock_start;
  timer::time_point clock_stop;
  std::chrono::duration<double> duration;

  // Start the chrono timer
  clock_start = timer::now();

  if (_max_iterations != 0) {
    // Run the coarse propagator
    runCoarseSolver();

    for (paraIndex k = 0; k < _max_iterations; k++) {
      // loop over each index of the precomputed values (in parallel)
      for (timeIndex n = 0; n < _coarse_solver->getNumTimeSteps() - 1; n++) {
	runFineSolver(n);
      }

      // Update the coarse solution
      update(k);
    }
  }
  else {
    runCoarseSolver();
    _global_output = _coarse_solver->getSolution();
  }

  clock_stop = timer::now();

  duration = clock_stop - clock_start;

  // set solve time
  _global_output->setSolveTime(duration.count());
}

template <typename Coarse, typename Fine>
void Parareal<Coarse, Fine>::writeToXML(pugi::xml_document& doc) const {
  std::ofstream out(_outpath);

  // Create teh root level node
  doc.append_child("parareal");

  // Write to xml
  _global_output->writeToXML(doc);

  // Save xml document
  doc.save(out);
}
