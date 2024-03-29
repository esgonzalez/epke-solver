#include <iostream>

#include "../catch.hpp"
#include "parareal/definitions.hpp"
#include "parareal/parareal.hpp"
#include "parareal/solver_parameters.hpp"
#include "parareal/solver_output.hpp"
#include "parareal/precursor.hpp"

TEST_CASE("Test parareal functions.", "[parareal]") {
  using namespace para;

  // EPKE parameters
  timeBins time = {0.0, 1.0, 2.0, 3.0};
  precBins<Precursor::ptr> precursors;
  timeBins rho_imp, gen_time, pow_norm, beta_eff, lambda_h;
  double theta, gamma_d, eta;

  // EPKE precomputed values
  timeBins precomp_time;
  timeBins power;
  timeBins rho;
  precBins<timeBins> concentrations;

  auto precomputed = std::make_shared<epke::EPKEOutput>(precomp_time,
							concentrations,
							power,
							rho);

  // additional parareal parameters
  precIndex max_iterations;
  std::string outpath;

  SECTION("Generate a fine time mesh for a given index", "[generateFineTime]") {
    EPKEParameters::ptr parameters =
      std::make_shared<EPKEParameters>(time,
				       precursors,
				       rho_imp,
				       gen_time,
				       pow_norm,
				       beta_eff,
				       lambda_h,
				       theta,
				       gamma_d,
				       eta);

    // build EPKE solver
    epke::Propagator solver(parameters, precomputed);

    timeIndex n_fine_per_coarse = 2;

    // build parareal solver
    Parareal parareal(solver,
		      outpath,
		      max_iterations,
		      n_fine_per_coarse);

    // test time interpolation
    timeBins fine_time_0 = {0.0, 0.5, 1.0};
    timeBins fine_time_1 = {0.0, 0.5, 1.0, 1.5, 2.0};
    timeBins fine_time_2 = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

    REQUIRE(parareal.generateFineTime(0) == fine_time_0);
    REQUIRE(parareal.generateFineTime(1) == fine_time_1);
    REQUIRE(parareal.generateFineTime(2) == fine_time_2);
  }

  SECTION("Generate fine time with pre-interpolated parameters",
	  "[generateFineTime]") {
    EPKEParameters::ptr parameters =
      std::make_shared<EPKEParameters>(time,
				       precursors,
				       rho_imp,
				       gen_time,
				       pow_norm,
				       beta_eff,
				       lambda_h,
				       theta,
				       gamma_d,
				       eta,
				       true);

    // build EPKE solver
    epke::Propagator solver(parameters, precomputed);

    timeIndex n_fine_per_coarse = 2;

    // build parareal solver
    Parareal parareal(solver,
		      outpath,
		      max_iterations,
		      n_fine_per_coarse);

    REQUIRE(parareal.generateFineTime(0) == time);
    REQUIRE(parareal.generateFineTime(1) == time);
    REQUIRE(parareal.generateFineTime(2) == time);
  }
}
