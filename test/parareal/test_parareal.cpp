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

    timeBins coarse_time = {0.0, 1.0, 2.0, 3.0};
    precBins<Precursor::ptr> coarse_precursors;

    SolverParameters coarse_parameters(coarse_time, coarse_precursors);
    SolverOutput coarse_output;
    precIndex K = 5;
    timeIndex n_fine_per_coarse = 2;

    Parareal parareal(coarse_parameters,
		      coarse_output,
		      K,
		      n_fine_per_coarse);

    timeBins fine_time_0 = {0.0, 0.5, 1.0};
    timeBins fine_time_1 = {0.0, 1.0, 1.5, 2.0};
    timeBins fine_time_2 = {0.0, 1.0, 2.0, 2.5, 3.0};

    REQUIRE(parareal.generateFineTime(0) == fine_time_0);
    REQUIRE(parareal.generateFineTime(1) == fine_time_1);
    REQUIRE(parareal.generateFineTime(2) == fine_time_2);
  }
}
