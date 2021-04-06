#ifndef _PARAREAL_SOLVER_PARAMETERS_HEADER_
#define _PARAREAL_SOLVER_PARAMETERS_HEADER_

class SolverParameters {
protected:
  template <typename T>
  using precBins = std::vector<T>;      // binning over precursor groups
  using timeBins = std::vector<double>; // binning over time variable
  using timeIndex = uint32_t;
  using precIndex = uint8_t;
};

class EPKEParameters : SolverParameters {
private:
  // coarse time mesh parameters
  const timeBins _gen_time; // Mean neutron generation time (Lambda)
  const timeBins _lambda_h; // Linear heat conduction constant
  const timeBins _pow_norm; // Power normalization factor
  const precBins<Precursor::ptr> _precursors;
  
public:
  EPKEParameters(const timeBins& gen_time,
		 const timeBins& lambda_h,
		 const timeBins& pow_norm,
		 const precBins<Precursor::ptr> precursors)
    : SolverParameters(),
      _gen_time(gen_time),
      _lambda_h(lambda_h),
      _pow_norm(pow_norm),
      _precursors(precursors) {}
};

#endif
