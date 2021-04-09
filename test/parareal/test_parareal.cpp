#include <iostream>

#include "../catch.hpp"
#include "parareal/parareal.hpp"

TEST_CASE("Test parareal functions.", "[parareal]") {
  SECTION("Generate a fine time mesh for a given index", "[generateFineTime]") {
    Parareal::TimeBins coarse_time = {0.0, 1.0, 2.0, 3.0};
    Parareal::TimeBins coarse_solution(coarse_time.size());
    int K = 5;
    int n_fine_per_coarse = 2;

    Parareal parareal(coarse_time, coarse_solution, K, n_fine_per_coarse);

    Parareal::TimeBins fine_time_0 = {0.0, 0.5, 1.0};
    Parareal::TimeBins fine_time_1 = {0.0, 1.0, 1.5, 2.0};
    Parareal::TimeBins fine_time_2 = {0.0, 1.0, 2.0, 2.5, 3.0};

    REQUIRE(parareal.generateFineTime(0) == fine_time_0);
    REQUIRE(parareal.generateFineTime(1) == fine_time_1);
    REQUIRE(parareal.generateFineTime(2) == fine_time_2);
  }
}
