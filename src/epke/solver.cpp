#include "epke/solver.hpp"
#include "utility/interpolate.hpp"

using namespace epke;

Solver::Solver(const Params::ptr params, const Output::ptr solution)
  : para::Solver(params, solution), _params(params), _solution(solution) {}

const double Solver::computeDT(const timeIndex n) const {
  return n > 0 ? _params->getTime(n) - _params->getTime(n-1) :
    _params->getTime(n+1) - _params->getTime(n);
}

const double Solver::computeGamma(const timeIndex n) const {
  return n < 2 ? 1.0 : computeDT(n-1) / computeDT(n);
}

const double Solver::computeOmega(const precIndex j,
				  const timeIndex n) const {
  using util::omegaN;

  const auto lj    = _params->getDecayConstant(j,n);
  const auto dt    = computeDT(n);
  const auto gamma = computeGamma(n);

  return _params->getGenTime(0) / _params->getGenTime(n) *
    _params->getDelayedFraction(j,n) / lj * omegaN(lj, dt, gamma);
}

const double Solver::computeZetaHat(const precIndex j,
				    const timeIndex n) const {
  using util::E, util::omegaN1, util::omegaN2;

  double beta_prev_prev, power_prev_prev, gen_time_prev_prev;

  const auto lj    = _params->getDecayConstant(j,n);
  const auto dt    = computeDT(n);
  const auto gamma = computeGamma(n);

  if (n < 2) {
    beta_prev_prev = 0.;
    power_prev_prev = 0.;
    gen_time_prev_prev = _params->getGenTime(0);
  } else {
    beta_prev_prev = _params->getDelayedFraction(j,n-2);
    power_prev_prev = getPower(n-2);
    gen_time_prev_prev = _params->getGenTime(n-2);
  }

  return E(lj, dt) * getConcentration(j, n-1) +
    1. / lj * _params->getGenTime(0) * getPower(n - 1) *
    _params->getDelayedFraction(j,n-1) / _params->getGenTime(n - 1) *
    omegaN1(lj, dt, gamma) +
    1. / lj * _params->getGenTime(0) * power_prev_prev * beta_prev_prev /
    gen_time_prev_prev * omegaN2(lj, dt, gamma);
}

const double Solver::computeA1(const timeIndex n) const {
  using util::omegaN;

  // define local variables to avoid numerous function calls
  const auto lh    = _params->getLambdaH(n);
  const auto dt    = computeDT(n);
  const auto gamma = computeGamma(n);

  return _params->getGammaD() * _params->getPowNorm(n) / lh *
    omegaN(lh, dt, gamma);
}

const double Solver::computeB1(const timeIndex n) const {
  using util::k0, util::E, util::omegaN1, util::omegaN2;

  // define local variables to avoid numerous function calls
  const auto lh    = _params->getLambdaH(n);
  const auto dt    = computeDT(n);
  const auto gamma = computeGamma(n);

  // we have to set these values so we don't get an out of range vector
  const auto H_prev_prev = n < 2 ? 0. : _params->getPowNorm(n-2) * getPower(n-2);

  return _params->getRhoImp(n) + E(lh, dt) * (getRho(n-1) -
					      _params->getRhoImp(n-1)) -
    1. / lh * getPower(0) * _params->getGammaD() *
    _params->getEta() * k0(lh, dt) +
    _params->getGammaD() / lh * (_params->getPowNorm(n-1) * getPower(n - 1) *
				 omegaN1(lh, dt, gamma) +
				 H_prev_prev * omegaN2(lh, dt, gamma));
}

const double Solver::computePower(const timeIndex n, const double alpha) const {
  const auto dt = computeDT(n);

  // accumulate the weighted sums
  double tau = 0.0, s_hat_d = 0.0, s_d_prev = 0.0;
  for (precIndex j = 0; j < _params->getNumPrecursors(); j++) {
    tau += _params->getDecayConstant(j,n) * computeOmega(j, n);
    s_hat_d += _params->getDecayConstant(j,n) * computeZetaHat(j, n);
    s_d_prev += _params->getDecayConstant(j,n-1) * getConcentration(j,n-1);
  }

  // compute the quadratic formula coefficients
  double a = _params->getTheta() * dt * computeA1(n) / _params->getGenTime(n);
  double b = _params->getTheta() * dt * (((computeB1(n) - _params->getBetaEff(n))
					/ _params->getGenTime(n) - alpha) +
				       tau / _params->getGenTime(0)) - 1;
  double c = _params->getTheta() * dt / _params->getGenTime(0) * s_hat_d +
    exp(alpha * dt) * ((1 - _params->getTheta()) * dt *
		       (((getRho(n - 1) - _params->getBetaEff(n - 1)) /
			 _params->getGenTime(n - 1) -
			 alpha) * getPower(n - 1) + s_d_prev /
			_params->getGenTime(0)) +
		       getPower(n - 1));

  // TODO: Add a quadratic formula util method to take care of this
  if (a < 0) {
    return (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
  } else if (a == 0) {
    return -c / b;
  } else {
    throw;
  }
}

void Solver::step(const timeIndex n) {
  // compute the transformation parameter
  double alpha = n > 1 ?
    1 / computeDT(n - 1) * log(getPower(n - 1) / getPower(n - 2)) : 0.0;

  // set solution time at this time step
  _solution->setTime(n, _params->getTime(n));

  // evaluate the power and normalization factor at this time step
  _solution->setPower(n, computePower(n, alpha));
  _solution->setPowNorm(n, _params->getPowNorm(n));

  // update the precursor concentrations
  for (precIndex j = 0; j < _params->getNumPrecursors(); j++) {
    _solution->setConcentration(j, n,
				getPower(n) * computeOmega(j,n) +
				computeZetaHat(j,n));
  }

  // update the full power and reactivity vectors
  _solution->setRho(n, computeA1(n) * getPower(n) + computeB1(n));
}

void Solver::solve() {
  // set initial conditions for the power and reactivity vectors
  for (timeIndex n = getStartTimeIndex(); n < getStopTimeIndex(); n++) {
    step(n);
  }
}
