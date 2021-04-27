#include <cmath>
#include <iomanip>
#include <sstream>

#include "epke_solver.hpp"
#include "utility/interpolate.hpp"

using namespace epke;

Solver::Solver(const pugi::xml_node& input_node,
	       const pugi::xml_node& output_node)
  : params(input_node),
    precomp(output_node),
    power(params.getNumTimeSteps()),
    rho(params.getNumTimeSteps()),
    concentrations(params.getNumPrecursors(),
		   timeBins(params.getNumTimeSteps())) {
  // set the initial concentration histories
  for (int k = 0; k < params.getNumPrecursors(); k++) {
    for (int n = 0; n < precomp.getNumTimeSteps(); n++) {
      concentrations[k][n] = precomp.getConcentration(k,n);
    }
  }

  // set the initial power histories and rho histories
  for (int n = 0; n < precomp.getNumTimeSteps(); n++) {
    power[n] = precomp.getPower(n);
    rho[n] = precomp.getRho(n);
  }
}

Solver::Solver(const EPKEParameters& params, const EPKEOutput& precomp)
  : params(params),
    precomp(precomp),
    power(params.getNumTimeSteps()),
    rho(params.getNumTimeSteps()),
    concentrations(params.getNumPrecursors(),
		   timeBins(params.getNumTimeSteps())) {
  // set the initial concentration histories
  for (int k = 0; k < params.getNumPrecursors(); k++) {
    for (int n = 0; n < precomp.getNumTimeSteps(); n++) {
      concentrations[k][n] = precomp.getConcentration(k,n);
    }
  }

  // set the initial power histories and rho histories
  for (int n = 0; n < precomp.getNumTimeSteps(); n++) {
    power[n] = precomp.getPower(n);
    rho[n] = precomp.getRho(n);
  }
}

Solver::ptr Solver::createFineSolver(const timeBins& fine_time,
				     const timeIndex coarse_index) {
  auto fine_solver =
    std::make_shared<Solver>(params.interpolate(fine_time),
			     precomp.createPrecomputed(coarse_index));
  // Push back the coarse solver's vector of fine solvers
  _fine_solvers.push_back(fine_solver);

  return fine_solver;
}

const para::SolverOutput::ptr Solver::assembleGlobalOutput() const {
  // TEMP
  return std::make_shared<para::SolverOutput>(precomp);
}

const double Solver::computeDT(const timeIndex n) const {
  return params.getTime(n) - params.getTime(n-1);
}

const double Solver::computeGamma(const timeIndex n) const {
  return n < 2 ? 1.0 : computeDT(n-1) / computeDT(n);
}

const double Solver::computeOmega(const precIndex k,
				  const timeIndex n) const {
  using util::k1, util::k2, util::E;

  const auto lk    = params.getDecayConstant(k,n);
  const auto dt    = computeDT(n);
  const auto gamma = computeGamma(n);
  const auto w     = 1. / E(lk, dt);

  return params.getGenTime(0) / params.getGenTime(n) * w *
    params.getDelayedFraction(k,n) * (k2(lk, dt) + gamma * dt *
				      k1(lk, dt)) / ((1 + gamma) * dt * dt);
}

const double Solver::computeZetaHat(const precIndex k,
				    const timeIndex n) const {
  using util::k0, util::k1, util::k2, util::E;
  double beta_prev_prev, power_prev_prev, gen_time_prev_prev;

  const auto lk    = params.getDecayConstant(k,n);
  const auto dt    = computeDT(n);
  const auto gamma = computeGamma(n);
  const auto w     = 1. / E(lk, dt);

  if (n < 2) {
    beta_prev_prev = params.getDelayedFraction(k,n-1);
    power_prev_prev = power.at(n-1);
    gen_time_prev_prev = params.getGenTime(n-1);
  } else {
    beta_prev_prev = params.getDelayedFraction(k,n-2);
    power_prev_prev = power.at(n-2);
    gen_time_prev_prev = params.getGenTime(n-2);
  }

  return w * concentrations.at(k).at(n - 1) +
    w * params.getGenTime(0) * power.at(n - 1) *
    params.getDelayedFraction(k,n-1) / params.getGenTime(n - 1) *
    (k0(lk, dt) - (k2(lk, dt) - dt * (gamma - 1) * k1(lk, dt)) /
     (gamma * dt * dt)) +
    w * params.getGenTime(0) * power_prev_prev * beta_prev_prev /
    gen_time_prev_prev * (k2(lk, dt) - dt * k1(lk, dt)) /
    ((1 + gamma) * gamma * dt * dt);
}

const double Solver::computeA1(const timeIndex n) const {
  using util::k1, util::k2, util::E;

  // define local variables to avoid numerous function calls
  const auto lh    = params.getLambdaH(n);
  const auto dt    = computeDT(n);
  const auto gamma = computeGamma(n);

  return params.getGammaD() * params.getPowNorm(n) / E(lh, dt) *
    (k2(lh, dt) + k1(lh, dt) * gamma * dt) / ((1 + gamma) * dt * dt);
}

const double Solver::computeB1(const timeIndex n) const {
  using util::k0, util::k1, util::k2, util::E;

  // define local variables to avoid numerous function calls
  const auto lh    = params.getLambdaH(n);
  const auto dt    = computeDT(n);
  const auto gamma = computeGamma(n);

  // we have to set these values so we don't get an out of range vector
  const auto n_accum = n < 2 ? n - 1 : n - 2;
  const auto H_prev_prev = params.getPowNorm(n_accum) * power.at(n_accum);

  return params.getRhoImp(n) + 1 / E(lh, dt) *
    ((rho.at(n-1) - params.getRhoImp(n-1)) - power.at(0) * params.getGammaD() *
     params.getEta() * k0(lh, dt)) + params.getGammaD() / E(lh, dt) *
    (params.getPowNorm(n-1) * power.at(n - 1) *
     (k0(lh, dt) - (k2(lh, dt) + (gamma - 1) * dt * k1(lh, dt)) /
      (gamma * dt * dt)) + H_prev_prev * (k2(lh, dt) - k1(lh, dt) * dt) /
     ((1 + gamma) * gamma * dt * dt));
}

const double Solver::computePower(const timeIndex n, const double alpha) const {
  const auto dt = computeDT(n);

  double tau = 0.0, s_hat_d = 0.0, s_d_prev = 0.0;
  for (int k = 0; k < params.getNumPrecursors(); k++) {
    // accumulate the weighted sum
    tau += params.getDecayConstant(k,n) * computeOmega(k, n);
    s_hat_d += params.getDecayConstant(k,n) * computeZetaHat(k, n);
    s_d_prev += params.getDecayConstant(k,n-1) * concentrations.at(k).at(n-1);
  }

  double a = params.getTheta() * dt * computeA1(n) / params.getGenTime(n);
  double b = params.getTheta() * dt * (((computeB1(n) - params.getBetaEff(n))
					/ params.getGenTime(n) - alpha) +
				       tau / params.getGenTime(0)) - 1;
  double c = params.getTheta() * dt / params.getGenTime(0) * s_hat_d +
    exp(alpha * dt) * ((1 - params.getTheta()) * dt *
		       (((rho.at(n - 1) - params.getBetaEff(n - 1)) /
			 params.getGenTime(n - 1) -
			 alpha) * power.at(n - 1) + s_d_prev /
			params.getGenTime(0)) +
		       power.at(n - 1));

  if (a < 0) {
    return (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
  } else if (a == 0) {
    return -c / b;
  } else {
    throw;
  }
}

// returns true if the transformation criterion is met
const bool Solver::acceptTransformation(const timeIndex n,
					const double alpha) const {

  const double power_prev_prev = n < 2 ? power.at(n - 1) : power.at(n - 2);

  double lhs = fabs(power.at(n) - exp(alpha * computeDT(n)) * power.at(n - 1));
  double rhs = fabs(power.at(n) - power.at(n - 1) -
		    (power.at(n - 1) - power_prev_prev) / computeGamma(n));
  return (lhs <= rhs);
}

para::SolverOutput::ptr Solver::solve() {
  // set initial conditions for the power and reactivity vectors
  double alpha = 0.0;

  for (int n = precomp.getNumTimeSteps(); n < params.getNumTimeSteps(); n++) {
    // compute the transformation parameter
    if (n > 1) {
      alpha = 1 / computeDT(n - 1) * log(power.at(n - 1) / power.at(n - 2));
    }

    // evaluate the power at this time step
    power[n] = computePower(n, alpha);

    // test whether we accept or reject the transformation parameter
    // if (!acceptTransformation(n, alpha, gamma)) {
    // alpha = 0.0;
    // power[n] = computeABC(n, alpha);
    //}

    // update the precursor concentrations
    for (int k = 0; k < params.getNumPrecursors(); k++) {
      concentrations[k][n] = power.at(n) * computeOmega(k,n) +
	computeZetaHat(k,n);
    }

    // update the full power and reactivity vectors
    rho[n] = computeA1(n) * power.at(n) + computeB1(n);
  }

  return std::make_shared<EPKEOutput>(power, rho, concentrations);
}

void Solver::buildXMLDoc(pugi::xml_document& doc) const {
  pugi::xml_node output_node = doc.append_child("epke_output");
  pugi::xml_node time_node = output_node.append_child("time");
  pugi::xml_node power_node = output_node.append_child("power");
  pugi::xml_node rho_node = output_node.append_child("rho");
  pugi::xml_node concs_node = output_node.append_child("concentrations");
  std::ostringstream time_str, power_str, rho_str, conc_str;

  for (int n = 0; n < params.getNumTimeSteps(); n++) {
    time_str << std::setprecision(6) << params.getTime(n);
    power_str << std::setprecision(12) << params.getPowNorm(n) * power.at(n);
    rho_str << std::setprecision(12) << rho.at(n);
    if (n != params.getNumTimeSteps() - 1) {
      time_str << " ";
      power_str << " ";
      rho_str << " ";
    }
  }

  time_node.text() = time_str.str().c_str();
  power_node.text() = power_str.str().c_str();
  rho_node.text() = rho_str.str().c_str();

  for (int k = 0; k < params.getNumPrecursors(); k++) {
    conc_str.str("");
    conc_str.clear();
    pugi::xml_node conc_node = concs_node.append_child("concentration");
    conc_node.append_attribute("k") = k;
    for (int n = 0; n < params.getNumTimeSteps(); n++) {
      conc_str << std::setprecision(12) << concentrations.at(k).at(n);

      if (n != params.getNumTimeSteps() - 1) {
        conc_str << " ";
      }
    }
    conc_node.text() = conc_str.str().c_str();
  }
}
