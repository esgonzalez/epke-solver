#ifndef _PARAREAL_PRECURSOR_HEADER_
#define _PARAREAL_PRECURSOR_HEADER_

#include <memory>
#include <vector>

template <typename T>
using timeBins = std::vector<T>;
using timeIndex = uint32_t;

class Precursor {
private:
  const timeBins<double> decay_constant; // lambda
  const timeBins<double> delayed_fraction; // beta
public:
  typedef std::shared_ptr<Precursor> ptr;
  Precursor(
    const timeBins<double>& decay_constant,
    const timeBins<double>& delayed_fraction )
    : decay_constant(decay_constant), delayed_fraction(delayed_fraction) {}

  // accessors
  const timeBins<double> &decayConstant() const { return decay_constant; }
  const timeBins<double> &delayedFraction() const { return delayed_fraction; }

  const double decayConstant(timeIndex n) const {
    return decay_constant.at(n);
  }

  const double delayedFraction(timeIndex n) const {
    return delayed_fraction.at(n);
  }
};

#endif
