#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <array>
#include <cmath>
#include <functional>
#include <limits>

extern "C" {

  // forward declaration of sgemm
  void
  dgesv_ (const int *N, const int *NRHS, double *A,
          const int *LDA, int* IPIV, double *B,
          const int *LDB, int *INFO);

}

using ratefunwjac = std::function<void (std::vector<double>::const_iterator starty0,
				    std::vector<double>::const_iterator endy0,
				    std::vector<double>::iterator startdydx,
				    std::vector<double> &jacobian,
				    bool compute_jacobian)>;

struct
rosenbrock_trapz {

  /*
    References: 

    [1] Hairer, Ernst; Nørsett, Syvert Paul; Wanner, Gerhard (2008),
    Solving ordinary differential equations I: Nonstiff problems,
    Berlin, New York: Springer-Verlag, ISBN 978-3-540-56670-0

    [2] Hairer, Ernst; Wanner, Gerhard,
    Solving ordinary differential equations II: Stiff and Differential-Algebraic problems,
    Berlin, New York: Springer-Verlag, ISBN 978-3-642-05220-0

  */

  std::vector<double> y_lo;
  std::vector<double> y_hi;
  std::vector<double> K0;
  std::vector<double> K1;
  std::vector<double> ys;
  
  double err;
  
  void
  step (std::vector<double>::const_iterator starty0,
	std::vector<double>::const_iterator endy0,
	double dx, ratefunwjac fun) {

    int ncomp = endy0 - starty0;
    static int ii = 0, jj = 0;
    std::vector<double> jacobian;

    jacobian.assign (ncomp*ncomp, 0.);
    y_lo.assign (ncomp, 0.);
    y_hi.assign (ncomp, 0.);
    K0.assign (ncomp, 0.);
    K1.assign (ncomp, 0.);
    ys.assign (ncomp, 0.);
     
    fun (starty0, endy0, K0.begin (), jacobian, false);
    
    std::copy (starty0, endy0, ys.begin ());
    auto itk = K0.begin ();
    for (auto & refys : ys)
      refys += *(itk++) * dx/2.;

    fun (ys.begin (), ys.end (), K1.begin (), jacobian, true);
    for (ii = 0; ii < ncomp; ++ii) {
      for (jj = 0; jj < ncomp; ++jj) {
	*(jacobian.begin () + ii + ncomp * jj) *= (- dx/2.);
	if (ii == jj)
	  *(jacobian.begin () + ii + ncomp * jj) += 1.;
      }
    }

    const int one = 1;
    int info = 0;
    std::vector<int> ipiv (ncomp, 0);
    dgesv_ (&ncomp, &one, jacobian.data (),
	    &ncomp, ipiv.data (), K1.data (),
	    &ncomp, &info);
    
    std::copy (starty0, endy0, y_lo.begin ());
    auto itk0 = K0.begin ();
    for (auto & refylo : y_lo)
      refylo += *(itk0++) * dx;

    std::copy (starty0, endy0, y_hi.begin ());
    itk0 = K0.begin ();
    auto itk1 = K1.begin ();
    for (auto & refyhi : y_hi)
      refyhi += (*(itk0++) + *(itk1++))* dx/2.;

    err = 0.;
    for (auto it_lo = y_lo.begin (), it_hi = y_hi.begin ();
	 it_lo != y_lo.end ();
	 ++it_lo, ++it_hi) {
      auto tmp = std::abs ((*it_lo) - (*it_hi));
      err = tmp > err ? tmp : err;
    }
  }

};


void
integrate_adaptive_trapz (ratefunwjac fun,
			  const std::vector<double> &y0,
			  double xstart, double xend, double tol,
			  std::vector<double> &x, std::vector<double> &y) {

  rosenbrock_trapz rk{};
  
  x.clear ();
  y.clear ();
  std::size_t ncomp = y0.end ()-y0.begin ();
    
  x.push_back (xstart);
  y.insert (y.end (), y0.begin (), y0.end ());
  auto y_back = y.begin ();
  
  double tnew = xstart;
  const double dtmin = 1.0e-14;
  double dt = 5*dtmin;
  
  int rejected = 0;
  while (x.back () < xend) {
    tnew = x.back () + dt;
    if (tnew > xend) {
      tnew = xend;
      dt = tnew - x.back ();
    }

    rk.step (y_back, y.end (), dt, fun);
    

    auto tmpfun  = [] (double x, double y) { return std::abs (x) < std::abs (y); };
    auto normylo = std::abs (*std::max_element (rk.y_lo.begin (), rk.y_lo.end (), tmpfun));

    std::cerr << "t = " << tnew << " dt = " << dt << " dtmin " << dtmin << std::endl;
    if (rk.err < tol * (1. + std::abs (normylo)) || dt <= dtmin) {
      x.push_back (tnew);
      y.insert (y.end (), rk.y_lo.begin (), rk.y_lo.end ());
      y_back = y.end () - ncomp;
     
    } else {
      ++rejected;
    }
    
    dt = .5 * dt * std::pow (tol / rk.err, 1./2.);
    
    if (dt < dtmin)
      dt = dtmin;
  }
  std::cerr << rejected << " rejected steps" << std::endl;
};

#endif
