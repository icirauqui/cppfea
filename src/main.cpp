#include <iostream>
#include <vector>
#include <math.h>

#include "element_test.cpp"
#include "simulation.cpp"



int main(int argc, char** argv) {
  //test_c3d6_displacements();
  //test_c3d6_loads();

  //test_c3d8_displacements();
  //test_c3d8_loads();
  
  //test_c2d4_displacements();
  //test_c2d4_loads();

  //simulation_manual();
  simulation_optimizer();



  
  return 0;
}
