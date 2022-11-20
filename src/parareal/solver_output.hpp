#ifndef _PARAREAL_SOLVER_OUTPUT_HEADER_
#define _PARAREAL_SOLVER_OUTPUT_HEADER_

#include <memory>

#include "parareal/definitions.hpp"

namespace pugi {
  class xml_document;
}

namespace para {

class SolverOutput {
public:
  using timeBins  = para::timeBins;
  using timeIndex = para::timeIndex;
  using ptr       = std::shared_ptr<SolverOutput>;

protected:
  // Time discretization
  timeBins _time;

  // Number of precomputed values (initial conditions)
  timeIndex _n_start;

  // Index to stop solve
  timeIndex _n_stop;

  double _solve_time;

public:
  // Construct empty from size
  SolverOutput(const timeIndex num_time_steps,
	       const timeIndex n_start,
	       const timeIndex n_stop)
    : _time(num_time_steps, 0.), _n_start(n_start), _n_stop(n_stop) {}

  // Construct from data vectors
  SolverOutput(const timeBins& time,
	       const timeIndex n_start,
	       const timeIndex n_stop)
    : _time(time), _n_start(n_start), _n_stop(n_stop) {}

  // Getters
  const timeIndex getNumTimeSteps() const { return _time.size(); }

  const timeIndex getStartTimeIndex() const { return _n_start; }

  const timeIndex getStopTimeIndex() const { return _n_stop; }

  const double getTime(const timeIndex n) const { return _time[n]; }

  // Set solve time
  void setSolveTime(double solve_time) { _solve_time = solve_time; }

  // Set the time at index n
  void setTime(const timeIndex n, const double val) { _time[n] = val; }

  // Create output object with truncated precomputed values from coarse solver
  virtual SolverOutput::ptr
  createPrecomputedImpl(const timeIndex n,
			const timeIndex n_fine_per_coarse) const = 0;

  // Create output object with only information from the coarse time steps
  //virtual SolverOutput::ptr coarsenImpl(const timeBins& coarse_time) const = 0;

  // Resize the number of time steps
  virtual void resize(const timeIndex n_steps) = 0;

  // Update the new coarse solution
  virtual void updateCoarseImpl(const timeIndex n, ptr coarse) = 0;

  // Update the solution with parareal adjustments
  virtual void updatePararealImpl(const timeIndex n,
				  ptr new_coarse,
				  ptr fine_coarsened,
				  ptr old_coarse) = 0;

  // Write the output to an xml document
  virtual void writeToXML(pugi::xml_document& doc) const = 0;
};

  template<typename T>
  std::shared_ptr<T> createPrecomputed(std::shared_ptr<T> precomp,
				       const timeIndex n,
				       const timeIndex n_fine_per_coarse) {
    return std::static_pointer_cast<T>(
	       precomp->createPrecomputedImpl(n,
					      n_fine_per_coarse));
  }

  template<typename T>
  std::shared_ptr<T> coarsen(std::shared_ptr<T> fine_output,
			     const timeBins& coarse_time) {
    return std::static_pointer_cast<T>(fine_output->coarsenImpl(coarse_time));
  }

  template<typename T>
  void updateCoarse(const timeIndex n,
		    std::shared_ptr<T> solution,
		    std::shared_ptr<T> coarse) {
    std::static_pointer_cast<T>(solution)->updateCoarseImpl(n, coarse);
  }

  template<typename T>
  void updateParareal(const timeIndex n,
		      std::shared_ptr<T> solution,
		      std::shared_ptr<T> new_coarse,
		      std::shared_ptr<T> fine_coarsened,
		      std::shared_ptr<T> old_coarse) {
    std::static_pointer_cast<T>(solution)->updatePararealImpl(n,
							      new_coarse,
							      fine_coarsened,
							      old_coarse);
  }

} // namespace para

#endif
