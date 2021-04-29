#include "../catch.hpp"
#include "parareal/definitions.hpp"
#include "parareal/solver_output.hpp"

TEST_CASE("Test solver_output functions.", "[SolverOutput]") {
  SECTION("Create precomputed solver output values", "[createPrecomputed]") {
    using namespace para;

    // output data
    timeBins time  = {0.0, 1.0, 2.0, 3.0, 4.0};
    timeBins power = {1.0, 2.0, 3.0, 2.5, 2.0};
    timeBins rho   = {0.0, 1.0, 2.0, 1.0, 0.0};
    precBins<timeBins> concentrations = {{0.5, 1.0, 1.5, 2.0, 2.5}};

    auto output = std::make_shared<epke::EPKEOutput>(time,
						     concentrations,
						     power,
						     rho);

    timeBins fine_p_history = {1.0, 2.0, 3.0};
    timeBins fine_rho_history = {0.0, 1.0, 2.0};
    precBins<timeBins> fine_zeta_histories = {{0.5, 1.0, 1.5}};

    auto precomp = createPrecomputed(output, 2);

    // Check the precomputed "history" values
    for (timeIndex n = 0; n < precomp->getNumTimeSteps(); n++) {
      REQUIRE(precomp->getPower(n) == fine_p_history.at(n));
      REQUIRE(precomp->getRho(n) == fine_rho_history.at(n));

      for (precIndex k = 0; k < precomp->getNumPrecursors(); k++) {
	REQUIRE(precomp->getConcentration(k,n) ==
		fine_zeta_histories.at(k).at(n));
      }
    }
  }

  SECTION("Coarsen values from a fine time grid", "[coarsen]") {
    using namespace para;

    timeBins fine_time  = {0.0, 0.5, 1.0, 1.5, 2.0};
    timeBins fine_power = {1.0, 2.0, 3.0, 2.5, 2.0};
    timeBins fine_rho   = {0.0, 1.0, 2.0, 1.0, 0.0};
    precBins<timeBins> fine_concentrations = {{0.5, 1.0, 1.5, 2.0, 2.5}};

    auto fine_output = std::make_shared<epke::EPKEOutput>(fine_time,
							  fine_concentrations,
							  fine_power,
							  fine_rho);

    timeBins coarse_time  = {0.0, 1.0, 2.0};
    timeBins coarse_power = {1.0, 3.0, 2.0};
    timeBins coarse_rho   = {0.0, 2.0, 0.0};
    precBins<timeBins> coarse_concentrations = {{0.5, 1.5, 2.5}};

    auto coarse_output = coarsen(fine_output, coarse_time);

    for (timeIndex n = 0; n < coarse_output->getNumTimeSteps(); n++) {
      REQUIRE(coarse_output->getPower(n) == coarse_power.at(n));
      REQUIRE(coarse_output->getRho(n) == coarse_rho.at(n));

      for (precIndex k = 0; k < coarse_output->getNumPrecursors(); k++) {
	REQUIRE(coarse_output->getConcentration(k,n) ==
		coarse_concentrations.at(k).at(n));
      }
    }
  }
}
