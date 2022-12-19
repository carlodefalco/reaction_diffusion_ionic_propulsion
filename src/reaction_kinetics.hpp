#ifndef REACTION_KINETICS_HPP
#define REACTION_KINETICS_HPP

#include <array>
#include <iostream>
#include <map>
#include <vector>
#include <json.hpp>

struct
reaction {
  
  std::map<std::string, int> reactants;
  std::map<std::string, int> products;
  std::array<double, 2> rate_coeffs;

};

struct
reaction_list {

  std::vector<reaction> data;

  void
  read (const std::string &filename) {
    std::ifstream jfile (filename);
    if (! jfile) {
    std::cerr << "error reading data file" << std::endl;
    }
    nlohmann::json jsn;

    jfile >> jsn;
    jfile.close ();

    for (auto & ii : jsn) {
      data.push_back (reaction{ii[0], ii[1], ii[2]});
    }
  }

  void
  print (std::ostream &os) {

    for (auto & ii : data) {

      auto jj0 = ii.reactants.begin ();
      os << jj0->second << " " << jj0->first;
      for (++jj0; jj0 != ii.reactants.end (); ++jj0)
	os  << " + " << jj0->second << " " << jj0->first;
      os << " = ";
      
      auto jj1 = ii.products.begin ();
      os << jj1->second << " " << jj1->first;
      for (++jj1; jj1 != ii.products.end (); ++jj1)
	os  << " + " << jj1->second << " " << jj1->first;

      os << std::endl << "rate coeffs: ";
      for (auto & jj : ii.rate_coeffs)
	os <<  jj << " ";
      os << std::endl << std::endl;
    }
  }
  
};

#endif
