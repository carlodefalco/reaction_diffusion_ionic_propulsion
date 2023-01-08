#include <iomanip>
#include <iostream>
#include <vector>
#include <reaction_kinetics.hpp>
#include <runge_kutta.hpp>


// Simulate the reaction system (Robertson problem)
//  A     -> B
//  2 B   -> C + B
//  B + C -> A + C
// with fwd coefficients .04, 3.e7 and 1.e4.
//
// reactions and coefficients are read from the
// file robertson_autocatalysis.json.
//
// Initial conditions A = 1, B = 0, C = 0.

reaction_list rl{};

void
fun (double x0,
     std::vector<double>::const_iterator starty0,
     std::vector<double>::const_iterator endy0,
     std::vector<double>::iterator startdydx)
{
  
  rl.change_rate (starty0, endy0, startdydx);
  
};


int
main () {

  dorpri5 rk{};
  rl.read ("robertson_autocatalysis.json");
  rl.pretty_print (std::cerr);
  
  double x0   = 0.;
  double xend = 1000.;


  std::vector<double> y0{1., 0., 0.};
  double tol  = 1.0e-8;

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
