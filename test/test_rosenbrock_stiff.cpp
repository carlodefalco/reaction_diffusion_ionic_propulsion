#include <iomanip>
#include <iostream>
#include <vector>
#include <rosenbrock_trapz.hpp>


// Simulate the reaction system (Robertson problem)
//  A     -> B
//  2 B   -> C + B
//  B + C -> A + C
// with fwd coefficients .04, 3.e7 and 1.e4.
// Initial conditions A = 1, B = 0, C = 0.

void
fun (std::vector<double>::const_iterator starty0,
     std::vector<double>::const_iterator endy0,
     std::vector<double>::iterator startdydx,
     std::vector<double> &jacobian,
     bool computejacobian = false)
{
  
  static double R1 = 0.;
  static double R2 = 0.;
  static double R3 = 0.;
  
  const double& A = *(starty0);
  const double& B = *(starty0+1);
  const double& C = *(starty0+2);

  double &Adot = *(startdydx);
  double &Bdot = *(startdydx+1);
  double &Cdot = *(startdydx+2);

  R1 = .04 * A;
  R2 = 3.e7 * B * B;
  R3 = 1.e4 * B*C;

  Adot = - R1 + R3;
  Bdot =   R1 - 2*R2 + R2 - R3;
  Cdot =   R2 - R3 + R3;

  if (computejacobian) {
    jacobian[0 + 3 * 0] = -.04;   jacobian[0 + 3 * 1] = 1e4 * C;               jacobian[0 + 3 * 0] =  1e4 * B; 
    jacobian[1 + 3 * 0] =  .04;   jacobian[1 + 3 * 1] = -2*3e7 * B - 1e4 * C;  jacobian[1 + 3 * 0] = -1e4 * B;
    jacobian[2 + 3 * 0] =    0;   jacobian[2 + 3 * 1] =  2*3e7 * B;            jacobian[2 + 3 * 0] =  0; 
  }
};


int
main () {

  rosenbrock_trapz rk{};

  
  double x0   = 0.;
  double xend = 1000.;


  std::vector<double> y0{1., 0., 0.};
  double tol  = 1.0e-8;

  std::vector<double> x{}, y{};
    
  integrate_adaptive_trapz (fun, y0, x0, xend, tol, x, y);
    
  for (int is = 0; is < x.size (); ++is) {
    std::cout << std::setprecision (16) << x[is];
    for (int ic = 0; ic < y0.size (); ++ic) {
      std::cout << "  " << y[y0.size ()*is + ic];
    }
    std::cout << std::endl;
  }
    
  return 0;
    
}
