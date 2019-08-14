#ifndef _EPKE_UTILITY_HEADER_
#define _EPKE_UTILITY_HEADER_

#include <cmath>

double E( const double lambda, const double delta_t ) {
  return exp(lambda * delta_t);
}

double k0( const double lambda, const double delta_t ) {
  return (E(lambda, delta_t) - 1) / lambda;
}

double k1( const double lambda, const double delta_t ) {
  return (delta_t * E(lambda, delta_t) - k0(lambda, delta_t)) / lambda;
}

double k2( const double lambda, const double delta_t ) {
  return (delta_t * delta_t *  E(lambda, delta_t) - 2 * k1(lambda, delta_t))
         / lambda;
}

#endif