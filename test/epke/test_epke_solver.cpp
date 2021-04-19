#include <string>

#include "../catch.hpp"
#include "epke/epke_solver.hpp"
#include "pugi/pugixml.hpp"

TEST_CASE("Regression tests for the EPKE solver.", "[EPKESolver]") {
  SECTION("Control rod ejection with feedback", "[generateFineTime]") {
    std::string ifile_name = "regression/cr_ejection_feedback_in.xml";
    std::string ofile_name = "regression/cr_ejection_feedback_out.xml";
    
    pugi::xml_document ifile, ofile;
    pugi::xml_parse_result iresult = ifile.load_file(ifile_name.c_str());
    pugi::xml_parse_result oresult = ofile.load_file(ofile_name.c_str());
    pugi::xml_node inode = ifile.child("epke_input");
    pugi::xml_node onode = ofile.child("epke_output");

    EPKEParameters params(inode);
    

  }
}
