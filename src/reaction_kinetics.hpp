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

void
to_json (nlohmann::json& j, const reaction& r) {
  j = nlohmann::json{{"reactants", r.reactants},
		     {"products", r.products},
		     {"rate_coeffs", r.rate_coeffs}};
}

void
from_json (const nlohmann::json& j, reaction& r) {
  j.at("reactants").get_to(r.reactants);
  j.at("products").get_to(r.products);
  j.at("rate_coeffs").get_to(r.rate_coeffs);
}


struct
reaction_list {

  std::vector<reaction> data;
  std::map<std::string, int> species;

  void
  init_species () {
    std::map<std::string, int>{}.swap (species);
    for (auto const & ii : data) {
      for (auto const & jj : ii.reactants) {
	species[jj.first] = 0;
      }
      for (auto const & jj : ii.products) {
	species[jj.first] = 0;
      }
    }
    int idx{0};
    for (auto & ii : species) {
      ii.second = idx++;
    }
  }
  
  void
  read (const std::string &filename) {
    std::ifstream jfile (filename);
    if (! jfile) {
      std::cerr << "error reading data file" << std::endl;
    }
    nlohmann::json jsn;

    jfile >> jsn;
    jfile.close ();

    data = jsn.get<std::vector<reaction>> ();
    init_species ();
  }

  void
  pretty_print (std::ostream &os) {

    os << "index of species" << std::endl;
    for (auto const & ii : species) {
      os << ii.second << ") " << ii.first << std::endl;
    }
    os << std::endl;
    
    os << "list of reactions" << std::endl;
    for (auto const & ii : data) {

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

  void
  print (std::ostream &os) {

    nlohmann::json jsn = data;
    os << jsn.dump(2) << std::endl;
    
  }
  
};


#endif
