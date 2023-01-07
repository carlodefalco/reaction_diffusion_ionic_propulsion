#include <iomanip>
#include <iostream>
#include <vector>
#include <runge_kutta.hpp>


// Simulate the reaction
//  2 A + B <-> C
// with fwd coefficient 2 and bwd coefficient 1.
// Initial conditions A = 1, B = 2, C = 0.
// Equilibrium is expected for A = sqrt (C/(2*B))

void
fun (double x0,
     std::vector<double>::const_iterator starty0,
     std::vector<double>::const_iterator endy0,
     std::vector<double>::iterator startdydx)
{
  
  static double Rf = 0.;
  static double Rb = 0.;

  const double& A = *(starty0);
  const double& B = *(starty0+1);
  const double& C = *(starty0+2);

  double &Adot = *(startdydx);
  double &Bdot = *(startdydx+1);
  double &Cdot = *(startdydx+2);

  Rf = 2. * (std::pow (A, 2) * B);
  Rb = C;

  Adot = - 2* Rf + 2 * Rb;
  Bdot = - Rf + Rb;
  Cdot =   Rf - Rb;
  
};


int
main () {

  dorpri5 rk{};

  
  double x0   = 0.;
  double xend = 2.;


  std::vector<double> y0{1., 2., 0.};
  double tol  = 1.0e-6;

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
