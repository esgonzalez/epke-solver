#include "parareal/input.hpp"

int main( int argc , char *argv[] ) {
  // read the command line
  if ( argc > 1 ) {
    auto input = Input(argv[1]);
    input.execute();
  }
  return 0;
}
