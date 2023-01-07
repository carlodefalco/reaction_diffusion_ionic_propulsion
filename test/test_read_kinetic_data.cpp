#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <json.hpp>
#include <reaction_kinetics.hpp>

int
main () {

  reaction_list rl{};
  rl.read ("reaction_list.json");
  rl.pretty_print (std::cout);
  rl.print (std::cout);
  
  return 0;
    
}
