#ifndef _PARAREAL_SOLVER_PARAMETERS_HEADER_
#define _PARAREAL_SOLVER_PARAMETERS_HEADER_

#include <vector>
#include <memory>

#include "parareal/definitions.hpp"

namespace pugi {
  class xml_document;
}

namespace para {

class SolverParameters {
public:
  using timeBins  = para::timeBins;
  using timeIndex = para::timeIndex;
  using precIndex = para::precIndex;
  using ptr       = std::shared_ptr<SolverParameters>;

protected:
  const timeBins _time; // time points

public:
  SolverParameters(const timeBins& time) : _time(time) {}

  // Getters
  virtual const precIndex getNumPrecursors() const = 0;
  const timeIndex getNumTimeSteps()    const { return _time.size(); }
  const double    getTime(timeIndex n) const { return _time.at(n);  }
  const timeBins& getTime()            const { return _time;        }

  virtual SolverParameters::ptr
  interpolateImpl(const timeBins& fine_time) const = 0;

  virtual void writeToXML(pugi::xml_document& doc) const = 0;
};

  template<typename T>
  std::shared_ptr<T> interpolate(std::shared_ptr<T> params,
				 const timeBins& fine_time) {
    return std::static_pointer_cast<T>(params->interpolateImpl(fine_time));
  }

} // namespace para

#endif
