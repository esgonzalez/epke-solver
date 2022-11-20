#ifndef _EPKE_PRECURSOR_HEADER_
#define _EPKE_PRECURSOR_HEADER_

#include <memory>

#include "parareal/definitions.hpp"

namespace epke {

  class Precursor {
  public:
    using timeBins  = para::timeBins;
    using timeIndex = para::timeIndex;
    using ptr       = std::shared_ptr<Precursor>;

  private:
    const timeBins _decay_constant; // lambda
    const timeBins _delayed_fraction; // beta
  public:
    Precursor(const timeBins& decay_constant, const timeBins& delayed_fraction )
      : _decay_constant(decay_constant), _delayed_fraction(delayed_fraction) {}

    // accessors
    const timeBins &decayConstant() const { return _decay_constant; }
    const timeBins &delayedFraction() const { return _delayed_fraction; }

    const double decayConstant(timeIndex n) const {
      return _decay_constant.at(n);
    }

    const double delayedFraction(timeIndex n) const {
      return _delayed_fraction.at(n);
    }
  }; // end class Precursor
} // end namespace epke

#endif
