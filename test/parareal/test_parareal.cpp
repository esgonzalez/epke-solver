#include <iostream>

#include "../catch.hpp"
#include "parareal/definitions.hpp"
#include "parareal/parareal.hpp"
#include "parareal/solver_parameters.hpp"
#include "parareal/solver_output.hpp"
#include "parareal/precursor.hpp"

TEST_CASE("Test parareal functions.", "[parareal]") {
  SECTION("Generate a fine time mesh for a given index", "[generateFineTime]") {
    using namespace para;

    // build EPKE parameters
    timeBins time = {0.0, 1.0, 2.0, 3.0};
    precBins<Precursor::ptr> precursors;
    timeBins rho_imp    = {0.0, 1.0, 2.0, 1.0, 0.0};
    timeBins gen_time   = {1.0, 1.0, 1.0, 1.0, 1.0};
    timeBins pow_norm   = {1.0, 1.0, 1.0, 1.0, 1.0};
    timeBins beta_eff   = {0.1, 0.1, 0.1, 0.1, 0.1};
    timeBins lambda_h   = {2.0, 2.0, 2.0, 2.0, 2.0};
    double theta   = 0.0;
    double gamma_d = 0.0;
    double eta     = 1.0;

    EPKEParameters parameters(time,
			      precursors,
			      rho_imp,
			      gen_time,
			      pow_norm,
			      beta_eff,
			      lambda_h,
			      theta,
			      gamma_d,
			      eta);

    // build EPKE precomputed values
    timeBins power;
    timeBins rho;
    precBins<timeBins> concentrations;

    epke::EPKEOutput precomputed(power, rho, concentrations);

    // build EPKE solver
    epke::Solver solver(parameters, precomputed);

    // additional parareal parameters
    precIndex max_iterations;
    timeIndex n_fine_per_coarse = 2;
    std::string outpath;

    Parareal parareal(solver,
		      outpath,
		      max_iterations,
		      n_fine_per_coarse);

    timeBins fine_time_0 = {0.0, 0.5, 1.0};
    timeBins fine_time_1 = {0.0, 1.0, 1.5, 2.0};
    timeBins fine_time_2 = {0.0, 1.0, 2.0, 2.5, 3.0};

    REQUIRE(parareal.generateFineTime(0) == fine_time_0);
    REQUIRE(parareal.generateFineTime(1) == fine_time_1);
    REQUIRE(parareal.generateFineTime(2) == fine_time_2);
  }
}
