#include "../catch.hpp"
#include "parareal/definitions.hpp"
#include "parareal/solver_parameters.hpp"
#include "parareal/precursor.hpp"

TEST_CASE("Test solver_parameters functions.", "[SolverParameters]") {
  using namespace para;

  // input data
  timeBins time       = {0.0, 1.0, 2.0, 3.0, 4.0};
  timeBins lambda_k   = {0.2, 0.2, 0.2, 0.2, 0.2};
  timeBins beta_eff_k = {0.1, 0.1, 0.1, 0.1, 0.1};
  precBins<Precursor::ptr> precursors =
    {std::make_shared<Precursor>(beta_eff_k, lambda_k)};
  timeBins rho_imp    = {0.0, 1.0, 2.0, 1.0, 0.0};
  timeBins gen_time   = {1.0, 1.0, 1.0, 1.0, 1.0};
  timeBins pow_norm   = {1.0, 1.0, 1.0, 1.0, 1.0};
  timeBins beta_eff   = {0.1, 0.1, 0.1, 0.1, 0.1};
  timeBins lambda_h   = {2.0, 2.0, 2.0, 2.0, 2.0};
  double theta   = 0.0;
  double gamma_d = 0.0;
  double eta     = 1.0;

  timeBins fine_time       = {0.0, 1.0, 2.0, 2.5, 3.0};
  timeBins fine_lambda_k   = {0.2, 0.2, 0.2, 0.2, 0.2};
  timeBins fine_beta_eff_k = {0.1, 0.1, 0.1, 0.1, 0.1};
  precBins<Precursor::ptr> fine_precursors =
    {std::make_shared<Precursor>(fine_beta_eff_k, fine_lambda_k)};
  timeBins fine_rho_imp    = {0.0, 1.0, 2.0, 1.5, 1.0};
  timeBins fine_gen_time   = {1.0, 1.0, 1.0, 1.0, 1.0};
  timeBins fine_pow_norm   = {1.0, 1.0, 1.0, 1.0, 1.0};
  timeBins fine_beta_eff   = {0.1, 0.1, 0.1, 0.1, 0.1};
  timeBins fine_lambda_h   = {2.0, 2.0, 2.0, 2.0, 2.0};

  SECTION("Interpolate epke solver parameters", "[interpolate]") {
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

    auto fine_params = parameters.interpolate(fine_time);

    // Check to make sure we interpolated all of the parameters correctly
    for (timeIndex n = 0; n < fine_params.getNumTimeSteps(); n++) {
      REQUIRE(fine_params.getTime(n)    == fine_time.at(n));
      REQUIRE(fine_params.getRhoImp(n)  == fine_rho_imp.at(n));
      REQUIRE(fine_params.getGenTime(n) == fine_gen_time.at(n));
      REQUIRE(fine_params.getPowNorm(n) == fine_pow_norm.at(n));
      REQUIRE(fine_params.getBetaEff(n) == fine_beta_eff.at(n));
      REQUIRE(fine_params.getLambdaH(n) == fine_lambda_h.at(n));

      for (precIndex k = 0; k < fine_params.getNumPrecursors(); k++) {
	REQUIRE(fine_params.getDecayConstant(k,n) ==
		fine_precursors.at(k)->decayConstant(n));
	REQUIRE(fine_params.getDelayedFraction(k,n) ==
		fine_precursors.at(k)->delayedFraction(n));
      }
    }

    // Check the solver parameters that aren't dependent on time
    REQUIRE(fine_params.getTheta()  == theta);
    REQUIRE(fine_params.getGammaD() == gamma_d);
    REQUIRE(fine_params.getEta()    == eta);
  }

  SECTION("Parameters have already been interpolated", "[interpolate]") {
    EPKEParameters parameters(time,
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

    auto fine_params = parameters.interpolate(fine_time);

    // Check to make sure that we didn't create a new parameters object
    for (timeIndex n = 0; n < fine_params.getNumTimeSteps(); n++) {
      REQUIRE(fine_params.getTime(n)    == time.at(n));
      REQUIRE(fine_params.getRhoImp(n)  == rho_imp.at(n));
      REQUIRE(fine_params.getGenTime(n) == gen_time.at(n));
      REQUIRE(fine_params.getPowNorm(n) == pow_norm.at(n));
      REQUIRE(fine_params.getBetaEff(n) == beta_eff.at(n));
      REQUIRE(fine_params.getLambdaH(n) == lambda_h.at(n));

      for (precIndex k = 0; k < fine_params.getNumPrecursors(); k++) {
	REQUIRE(fine_params.getDecayConstant(k,n) ==
		precursors.at(k)->decayConstant(n));
	REQUIRE(fine_params.getDelayedFraction(k,n) ==
		precursors.at(k)->delayedFraction(n));
      }
    }

    // Check the solver parameters that aren't dependent on time
    REQUIRE(fine_params.getTheta()  == theta);
    REQUIRE(fine_params.getGammaD() == gamma_d);
    REQUIRE(fine_params.getEta()    == eta);
  }
}
