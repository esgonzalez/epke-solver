#include "../catch.hpp"
#include "utility/interpolate.hpp"

using namespace util;

TEST_CASE( "Test linear interpolation functions", "[interpolate]" ) {
  std::vector<double> x = {1.0, 1.5, 2.0};
  std::vector<double> y = {2.0, 3.0, 4.0};

  SECTION("Linearly interpolate a single value") {
    double x_val = 1.25;

    REQUIRE( interpolate(x,y,x_val) == 2.5 );
  }

  SECTION("Linearly interpolate a vector of values") {
    std::vector<double> x_new = {1.0, 1.25, 1.35, 1.5, 1.65, 1.95, 2.0};
    std::vector<double> y_new = {2.0, 2.5, 2.7, 3.0, 3.3, 3.9, 4.0};

    REQUIRE( interpolate(x,y,x_new) == y_new );
  }

  SECTION("Linearly interpolate a vector of length 1") {
    x = {0.0};
    y = {1e-6};

    std::vector<double> x_new = {0.0};
    std::vector<double> y_new = {1e-6};

    REQUIRE( interpolate(x,y,x_new) == y_new );
  }
}
