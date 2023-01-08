#include <iomanip>
#include <iostream>
#include <vector>
#include <reaction_kinetics.hpp>
#include <runge_kutta.hpp>


// Simulate the reaction system studied in
// N. Balcon, G. Hagelaar, J. Boeuf,
// Numerical model of an argon atmospheric pressure rf discharge,
// IEEE transactions on plasma science 36 (5) (2008) 2782– 2787.

reaction_list rl{};

void
fun (double x0,
     std::vector<double>::const_iterator starty0,
     std::vector<double>::const_iterator endy0,
     std::vector<double>::iterator startdydx)
{ rl.change_rate (starty0, endy0, startdydx); };


int
main () {

  dorpri5 rk{};
  rl.read ("balcon_et_al_argon_ionization.json");
  rl.pretty_print (std::cerr);
  
  double x0   = 0.;
  double xend = 1.0e-6;

  std::size_t nspecies = rl.species.size ();
  std::vector<double> y0(nspecies, 0.0);

  // Set initial state
  y0[rl.species.at ("Ar")] = 2.5e+19;
  y0[rl.species.at ("e")] = 1.0e+6;
  y0[rl.species.at ("Ar+")] = 1.0e+6;
  y0[rl.species.at ("Ar2+")] = 1.0e+3;
  y0[rl.species.at ("Ar*")] = 1.0e+10;
  //y0[rl.species.at ("Ar2*")] = 1.0e+9;
  
  double tol  = 1.0e-2;

  std::vector<double> x{}, y{};
    
  integrate_adaptive (fun, y0, x0, xend, tol, x, y);
    
  for (int is = 0; is < x.size (); ++is) {
    std::cout << std::setprecision (16) << x[is];
    for (int ic = 0; ic < y0.size (); ++ic) {
      std::cout << "  " << y[y0.size ()*is + ic];
    }
    std::cout << std::endl;
  }
    
  return 0;
    
}
