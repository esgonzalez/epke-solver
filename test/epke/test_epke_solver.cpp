#include <string>
#include <iostream>

#include "../catch.hpp"
#include "epke/epke_solver.hpp"
#include "parareal/parareal.hpp"
#include "parareal/solver_parameters.hpp"
#include "pugi/pugixml.hpp"
#include "utility/load_data.hpp"

// TODO: Move regression tests to the test_parareal.cpp file
TEST_CASE("Regression tests for the EPKE solver.", "[EPKESolver]") {
  SECTION("Control rod ejection with feedback", "[buildXMLDoc]") {
    std::string ifile = "test/epke/regression/cr_ejection_feedback_in.xml";
    std::string ofile = "test/epke/regression/cr_ejection_feedback_out.xml";

    pugi::xml_document idoc, odoc;
    pugi::xml_parse_result iresult = idoc.load_file(ifile.c_str());
    pugi::xml_parse_result oresult = odoc.load_file(ofile.c_str());

    if (!iresult) {
      std::cout << iresult.description() << std::endl;
      throw;
    }

    if (!oresult) {
      std::cout << oresult.description() << std::endl;
      throw;
    }

    pugi::xml_node parareal_node = idoc.child("parareal");

    para::Parareal parareal(parareal_node);

    auto n_steps =
      parareal_node.child("epke_input").attribute("n_steps").as_int();

    // Load the output power from the regression test
    para::timeBins opower =
      util::loadVectorData(odoc.child("epke_output").child("power"), n_steps);

    // Run the parareal solver
    parareal.solve();

    pugi::xml_document tdoc;

    // Build the test output xml document
    parareal.writeToXML(tdoc);

    para::timeBins tpower =
      util::loadVectorData(tdoc.child("epke_output").child("power"), n_steps);

    REQUIRE(tpower.size() == opower.size());

    for (para::timeIndex n = 0; n < tpower.size(); n++) {
      REQUIRE(tpower.at(n) == Approx(opower.at(n)));
    }
  }
}
