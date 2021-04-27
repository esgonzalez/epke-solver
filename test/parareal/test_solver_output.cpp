#include "../catch.hpp"
#include "parareal/definitions.hpp"
#include "parareal/solver_output.hpp"

TEST_CASE("Test solver_output functions.", "[SolverOutput]") {
  SECTION("Create precomputed solver output values", "[createPrecomputed]") {
    using namespace para;

    // output data
    timeBins power = {1.0, 2.0, 3.0, 2.5, 2.0};
    timeBins rho   = {0.0, 1.0, 2.0, 1.0, 0.0};
    precBins<timeBins> concentrations = {{0.5, 1.0, 1.5, 2.0, 2.5}};

    epke::EPKEOutput output(power, rho, concentrations);

    timeBins fine_p_history = {1.0, 2.0, 3.0};
    timeBins fine_rho_history = {0.0, 1.0, 2.0};
    precBins<timeBins> fine_zeta_histories = {{0.5, 1.0, 1.5}};

    auto precomp = output.createPrecomputed(2);

    // Check the precomputed "history" values
    for (timeIndex n = 0; n < precomp.getNumTimeSteps(); n++) {
      REQUIRE(precomp.getPower(n) == fine_p_history.at(n));
      REQUIRE(precomp.getRho(n) == fine_rho_history.at(n));

      for (precIndex k = 0; k < precomp.getNumPrecursors(); k++) {
	REQUIRE(precomp.getConcentration(k,n) ==
		fine_zeta_histories.at(k).at(n));
      }
    }
  }
}
