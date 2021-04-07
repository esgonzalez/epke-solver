#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "epke_solver.hpp"
#include "utility/interpolate.hpp"

Solver::Solver(
    const timeBins& time, const timeBins& gen_time, const timeBins& pow_norm,
    const timeBins& rho_imp, const timeBins& beta_eff, const timeBins& lambda_h,
    const precBins<Precursor::ptr>& precursors)
    : time(time), gen_time(gen_time), pow_norm(pow_norm), rho_imp(rho_imp),
      beta_eff(beta_eff), lambda_h(lambda_h),
      precursors(precursors), delta_t(time.size()), power(time.size()),
      rho(time.size()), rho_d(time.size()), omega(precursors.size()),
      zeta_hat(precursors.size()),
      concentrations(precursors.size(), timeBins(time.size())) {

  // compute the delta_t vector
  for (int i = 0; i < time.size() - 1; i++) {
    delta_t[i] = time.at(i + 1) - time.at(i);
  }
  delta_t[time.size() - 1] =
      time.at(time.size() - 1) - time.at(time.size() - 2);
}

void Solver::setPrecomputedValues(const timeBins& power_history,
				  const precBins<timeBins>& zeta_histories) {
  p_history = power_history;
  
  for (int k = 0; k < precursors.size(); k++) {
    for (int n = 0; n < p_history.size(); n++) {
      concentrations[k][n] = zeta_histories[k][n];
    }
  }

  for (int n = 0; n < p_history.size(); n++) {
    power[n] = p_history[n];
  }
}

const double Solver::computeOmega(
    const PrecIndex& k, const TimeIndex& n, const double& w,
    const double& gamma) const {
  const auto beta_eff_k = &precursors.at(k)->delayedFraction();
  const auto lambda_k = &precursors.at(k)->decayConstant();
  return gen_time.at(0) / gen_time.at(n) * beta_eff_k->at(n) * w *
         (k2(lambda_k->at(n), delta_t.at(n)) +
          gamma * delta_t.at(n) * k1(lambda_k->at(n), delta_t.at(n))) /
         ((1 + gamma) * delta_t.at(n) * delta_t.at(n));
}

const double Solver::computeZetaHat(
    const PrecIndex& k, const TimeIndex& n, const double& w,
    const double& gamma) const {
  const auto beta_eff_k = &precursors.at(k)->delayedFraction();
  const auto lambda_k = &precursors.at(k)->decayConstant();
  double beta_prev_prev, power_prev_prev, gen_time_prev_prev;

  if (n < 2) {
    beta_prev_prev = beta_eff_k->at(n - 1);
    power_prev_prev = power.at(1);
    gen_time_prev_prev = gen_time.at(n - 1);
  } else {
    beta_prev_prev = beta_eff_k->at(n - 2);
    power_prev_prev = power.at(n - 2);
    gen_time_prev_prev = gen_time.at(n - 2);
  }

  return w * concentrations.at(k).at(n - 1) +
         w * gen_time.at(0) * power.at(n - 1) * beta_eff_k->at(n - 1) /
             gen_time.at(n - 1) *
             (k0(lambda_k->at(n), delta_t.at(n)) -
              (k2(lambda_k->at(n), delta_t.at(n)) -
               delta_t.at(n) * (gamma - 1) *
                   k1(lambda_k->at(n), delta_t.at(n))) /
                  (gamma * delta_t.at(n) * delta_t.at(n))) +
         w * gen_time.at(0) * power_prev_prev * beta_prev_prev /
             gen_time_prev_prev *
             (k2(lambda_k->at(n), delta_t.at(n)) -
              delta_t.at(n) * k1(lambda_k->at(n), delta_t.at(n))) /
             ((1 + gamma) * gamma * delta_t.at(n) * delta_t.at(n));
}

std::pair<double, double> Solver::computeA1B1(
    const TimeIndex& n, const double gamma_d, const double eta,
    const double& gamma) {
  double H_prev_prev;

  // we have to set these values so we don't get an out of range vector
  if (n < 2) {
    H_prev_prev = pow_norm.at(n - 1) * power.at(n - 1);
  } else {
    H_prev_prev = pow_norm.at(n - 2) * power.at(n - 2);
  }

  // define local variables to avoid numerous function calls
  double lambda_h_n = lambda_h.at(n);
  double delta_t_n = delta_t.at(n);

  double a1 = gamma_d * pow_norm.at(n) / E(lambda_h_n, delta_t_n) *
              (k2(lambda_h_n, delta_t_n) +
               k1(lambda_h_n, delta_t_n) * gamma * delta_t_n) /
              ((1 + gamma) * delta_t_n * delta_t_n);
  double b1 = rho_imp.at(n) +
              1 / E(lambda_h_n, delta_t_n) *
                  (rho_d.at(n - 1) -
                   power.at(0) * gamma_d * eta * k0(lambda_h_n, delta_t_n)) +
              gamma_d / E(lambda_h_n, delta_t_n) *
                  (pow_norm.at(n - 1) * power.at(n - 1) *
                       (k0(lambda_h_n, delta_t_n) -
                        (k2(lambda_h_n, delta_t_n) +
                         (gamma - 1) * delta_t_n * k1(lambda_h_n, delta_t_n)) /
                            (gamma * delta_t_n * delta_t_n)) +
                   H_prev_prev *
                       (k2(lambda_h_n, delta_t_n) -
                        k1(lambda_h_n, delta_t_n) * delta_t_n) /
                       ((1 + gamma) * gamma * delta_t_n * delta_t_n));;

  return std::make_pair(a1, b1);
}

const double Solver::computePower(
    const TimeIndex& n, const double& theta, const double& gamma_d,
    const double& eta, const double& alpha, const double& gamma) {
  double tau = 0.0, s_hat_d = 0.0, s_d_prev = 0.0;
  for (int k = 0; k < precursors.size(); k++) {
    double w = 1 / E(precursors.at(k)->decayConstant().at(n), delta_t.at(n));
    const auto decay_constant = &precursors.at(k)->decayConstant();

    omega[k] = computeOmega(k, n, w, gamma);
    zeta_hat[k] = computeZetaHat(k, n, w, gamma);

    // accumulate the weighted sum
    tau += decay_constant->at(n) * omega.at(k);
    s_hat_d += decay_constant->at(n) * zeta_hat.at(k);
    s_d_prev += decay_constant->at(n - 1) * concentrations.at(k).at(n - 1);
  }

  std::pair<double, double> a1b1 = computeA1B1(n, gamma_d, eta, gamma);

  return computeABC(n, theta, alpha, a1b1, tau, s_hat_d, s_d_prev);
}

const double Solver::computeABC(
    const TimeIndex& n, const double theta, const double& alpha,
    const std::pair<double, double>& a1b1, const double& tau,
    const double& s_hat_d, const double& s_d_prev) const {
  double a = theta * delta_t.at(n) * a1b1.first / gen_time.at(n);
  double b = theta * delta_t.at(n) *
                 (((a1b1.second - beta_eff.at(n)) / gen_time.at(n) - alpha) +
                  tau / gen_time.at(0)) -
             1;
  double c =
      theta * delta_t.at(n) / gen_time.at(0) * s_hat_d +
      exp(alpha * delta_t.at(n)) *
          ((1 - theta) * delta_t.at(n) *
               (((rho.at(n - 1) - beta_eff.at(n - 1)) / gen_time.at(n - 1) -
                 alpha) *
                    power.at(n - 1) +
                s_d_prev / gen_time.at(0)) +
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
const bool Solver::acceptTransformation(
    const TimeIndex& n, const double& alpha, const double& gamma) const {

  const double power_prev_prev = n < 2 ? power.at(n - 1) : power.at(n - 2);

  double lhs = fabs(power.at(n) - exp(alpha * delta_t.at(n)) * power.at(n - 1));
  double rhs = fabs(
      power.at(n) - power.at(n - 1) -
      (power.at(n - 1) - power_prev_prev) / gamma);
  return (lhs <= rhs);
}

void Solver::solve(const double theta, const double gamma_d, const double eta) {
  std::cout << "solving..." << std::endl;

  // set initial conditions for the power and reactivity vectors
  rho[0] = rho_imp.at(0);
  double alpha = 0.0;

  // set the initial condition for the precursor concentrations
  for (int k = 0; k < precursors.size(); k++) {
    auto precursor = precursors.at(k);
    concentrations[k][0] = precursor->delayedFraction().at(0) * power.at(0) /
                           precursor->decayConstant().at(0);
  }

  for (int n = p_history.size(); n < time.size(); n++) {
    double gamma = delta_t.at(n - 1) / delta_t.at(n);

    // compute the transformation parameter
    if (n > 1) {
      alpha = 1 / delta_t.at(n - 1) * log(power.at(n - 1) / power.at(n - 2));
    }

    // evaluate the power at this time step
    power[n] = computePower(n, theta, gamma_d, eta, alpha, gamma);

    // test whether we accept or reject the transformation parameter
    // if (!acceptTransformation(n, alpha, gamma)) {
    // alpha = 0.0;
    // power[n] = computeABC(n, theta, alpha);
    //}

    // update the precursor concentrations
    for (int k = 0; k < precursors.size(); k++) {
      concentrations[k][n] = power.at(n) * omega.at(k) + zeta_hat.at(k);
    }

    std::pair<double, double> a1b1 = computeA1B1(n, gamma_d, eta, gamma);

    // update the full power and reactivity vectors
    rho[n] = a1b1.first * power.at(n) + a1b1.second;
    rho_d[n] = rho.at(n) - rho_imp.at(n);
  }

  std::cout << "completed solve" << std::endl;
}

void Solver::buildXMLDoc(pugi::xml_document& doc) const {
  pugi::xml_node output_node = doc.append_child("epke_output");
  pugi::xml_node time_node = output_node.append_child("time");
  pugi::xml_node power_node = output_node.append_child("power");
  pugi::xml_node rho_node = output_node.append_child("rho");
  pugi::xml_node concs_node = output_node.append_child("concentrations");
  std::ostringstream time_str, power_str, rho_str, conc_str;

  for (int n = 0; n < time.size(); n++) {
    time_str << std::setprecision(6) << time.at(n);
    power_str << std::setprecision(12) << pow_norm.at(n) * power.at(n);
    rho_str << std::setprecision(12) << rho.at(n);
    if (n != time.size() - 1) {
      time_str << ", ";
      power_str << ", ";
      rho_str << ", ";
    }
  }

  time_node.append_attribute("values") = time_str.str().c_str();
  power_node.append_attribute("values") = power_str.str().c_str();
  rho_node.append_attribute("values") = rho_str.str().c_str();

  for (int k = 0; k < precursors.size(); k++) {
    conc_str.str("");
    conc_str.clear();
    pugi::xml_node conc_node = concs_node.append_child("concentration");
    conc_node.append_attribute("k") = k;
    for (int n = 0; n < time.size(); n++) {
      conc_str << std::setprecision(12) << concentrations.at(k).at(n);

      if (n != time.size() - 1) {
        conc_str << ", ";
      }
    }
    conc_node.append_attribute("values") = conc_str.str().c_str();
  }
}
