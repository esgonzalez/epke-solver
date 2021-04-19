#ifndef _PARAREAL_PRECURSOR_HEADER_
#define _PARAREAL_PRECURSOR_HEADER_

#include <memory>

#include "parareal/definitions.hpp"
#include "pugi/pugixml.hpp"
#include "utility/load_data.hpp"

class Precursor {
public:
  using timeBins  = para::timeBins;
  using timeIndex = para::timeIndex;
  using ptr       = std::shared_ptr<Precursor>;

private:
  const timeBins _decay_constant; // lambda
  const timeBins _delayed_fraction; // beta
public:
  // constructor from pugi xml node
  Precursor(const pugi::xml_node& prec_node, timeIndex n_steps) :
    _decay_constant(util::loadVectorData(prec_node.child("decay_constant"), n_steps)),
    _delayed_fraction(util::loadVectorData(prec_node.child("delayed_fraction"), n_steps))
  {}

  // constructor from data vectors
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
};

#endif
