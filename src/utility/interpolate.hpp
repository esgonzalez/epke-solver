#ifndef _UTILITY_INTERPOLATE_HEADER_
#define _UTILITY_INTERPOLATE_HEADER_

#include <cmath>
#include <vector>

namespace util {

inline const double interpolate(const std::vector<double> &x,
				const std::vector<double> &y,
				const double x_val)
{
   int size = x.size();

   if (size == 1) {
     return y[0];
   }

   int i = 0;

   // find the index of the lower-bound x-value in the x-vector
   if (x_val >= x[size - 2]) {
     i = size - 2;
   }
   else {
     while (x_val > x[i+1]) { i++; }
   };

   // interpolation points
   double a = x[i], ya = y[i], b = x[i+1], yb = y[i+1];

   // slope
   double t = ( yb - ya ) / ( b - a );

   return ya + t * ( x_val - a );
}

inline const std::vector<double> interpolate(const std::vector<double> &x,
					     const std::vector<double> &y,
					     const std::vector<double> &x_new) {
  std::vector<double> y_new;

  for (const auto& x_val : x_new) {
    y_new.push_back(interpolate(x,y,x_val));
  }

  return y_new;
}

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int n_steps) {
  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(n_steps);

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end);

  return linspaced;

}

  inline double E(const double lambda, const double delta_t) {
    return exp(-lambda * delta_t);
  }

  inline double k0(const double lambda, const double delta_t) {
    return 1. - E(lambda, delta_t);
  }

  inline double k1(const double lambda, const double delta_t) {
    return 1. - k0(lambda, delta_t) / (lambda * delta_t);
  }

  inline double k2(const double lambda, const double delta_t) {
    return 1. - 2. * k1(lambda, delta_t) / (lambda * delta_t);
  }

  inline double omega0(const double lambda, const double delta_t) {
    return E(lambda, delta_t);
  }

  inline double omegaN(const double lambda,
		       const double delta_t,
		       const double gamma) {
    return (k2(lambda, delta_t) + gamma * k1(lambda, delta_t)) / (1 + gamma);
  }

  inline double omegaN1(const double lambda,
			const double delta_t,
			const double gamma) {
    return k0(lambda, delta_t) - (k2(lambda, delta_t) + (gamma - 1) *
				  k1(lambda, delta_t)) / gamma;
  }

  inline double omegaN2(const double lambda,
			const double delta_t,
			const double gamma) {
    return (k2(lambda, delta_t) - k1(lambda, delta_t)) / ((1 + gamma) * gamma);
  }

} // namespace util

#endif
