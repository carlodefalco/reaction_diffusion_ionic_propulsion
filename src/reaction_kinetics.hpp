#ifndef REACTION_KINETICS_HPP
#define REACTION_KINETICS_HPP

#include <array>
#include <fstream>
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


  /*

    given the set of M reactions involving N species

    \sum_{k=1}^{N} r_{i,k} s_k <--> \sum_{k=1}^{N} p_{i,k} s_k    i=1..M
    
    with forward reaction coefficient cf_i and backward reaction coefficient cb_i

    the rate Rf_i of the i-th forward reaction is 

    Rf_i = cf_i \prod_k s_k^{\left( r_{i,k} \right)}

    the rate Rb_i of the i-th forward reaction is 

    Rb_i = cb_i \prod_k s_k^{\left( p_{i,k} \right)}

    the net rate of the i-th reaction is 

    R_i = Rf_i - Rb_i

    the system of odes is

    \dot{s}_k = \sum_{i=1}^M \left( - Rf_i r_{i,k} + Rf_i p_{i,k}
                - Rb_i p_{i,k} - Rb_i r_{i,k} \right) =
                - \sum_{i=1}^M \left( R_i (r_{i,k} - p_{i,k}) \right) 

   */

  void
  change_rate (std::vector<double>::const_iterator state_begin,
	       std::vector<double>::const_iterator state_end,
	       std::vector<double>::iterator dot_state_begin) {

    // assuming the range
    // dot_state_begin ... dot_state_begin + (state_end-state_begin)
    // to be allocated and assigned with zeros
    
    auto dot_state_it = dot_state_begin;
    double Rfi = 0., Rbi = 0., Ri = 0.;
     
    for (auto reaction_it = data.begin ();
	 reaction_it != data.end ();
	 ++reaction_it, ++dot_state_it) {

      Rfi = reaction_it->rate_coeffs[0];
      for (auto const & ir : reaction_it->reactants) {
	int sidx = species.at (ir.first);
	Rfi *= std::pow (*(state_begin+sidx), ir.second);
      }

      Rbi = reaction_it->rate_coeffs[1];
      for (auto const & ir : reaction_it->products) {
	int sidx = species.at (ir.first);
	Rbi *= std::pow (*(state_begin+sidx), ir.second);
      }

      Ri = Rfi - Rbi;

      for (auto const & ir : reaction_it->reactants) {
	int sidx = species.at (ir.first);
	*(dot_state_begin+sidx) += - Ri * ir.second;
      }

       for (auto const & ir : reaction_it->products) {
	int sidx = species.at (ir.first);
	*(dot_state_begin+sidx) += Ri * ir.second;
      }
       
    }
	

  }


  void
  change_rate_wjac (std::vector<double>::const_iterator state_begin,
		    std::vector<double>::const_iterator state_end,
		    std::vector<double>::iterator dot_state_begin,
		    std::vector<double> &jacobian,
		    bool computejacobian) {

    // assuming the range
    // dot_state_begin ... dot_state_begin + (state_end-state_begin)
    // to be allocated and assigned with zeros
    
    auto dot_state_it = dot_state_begin;
    double Rfi = 0., Rbi = 0., Ri = 0.;
     
    for (auto reaction_it = data.begin ();
	 reaction_it != data.end ();
	 ++reaction_it, ++dot_state_it) {

      Rfi = reaction_it->rate_coeffs[0];
      for (auto const & ir : reaction_it->reactants) {
	int sidx = species.at (ir.first);
	Rfi *= std::pow (*(state_begin+sidx), ir.second);
      }

      Rbi = reaction_it->rate_coeffs[1];
      for (auto const & ir : reaction_it->products) {
	int sidx = species.at (ir.first);
	Rbi *= std::pow (*(state_begin+sidx), ir.second);
      }

      Ri = Rfi - Rbi;

      for (auto const & ir : reaction_it->reactants) {
	int sidx = species.at (ir.first);
	*(dot_state_begin+sidx) += - Ri * ir.second;
      }

       for (auto const & ir : reaction_it->products) {
	int sidx = species.at (ir.first);
	*(dot_state_begin+sidx) += Ri * ir.second;
      }
       
    }
	

  }

  
};


#endif
