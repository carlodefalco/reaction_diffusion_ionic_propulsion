#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <array>
#include <cmath>
#include <functional>
#include <limits>

using ratefun = std::function<void (double x0,
				    std::vector<double>::const_iterator starty0,
				    std::vector<double>::const_iterator endy0,
				    std::vector<double>::iterator startdydx)>;

struct
dorpri5 {

  /*
    References: 

    [1] Hairer, Ernst; Nørsett, Syvert Paul; Wanner, Gerhard (2008),
    Solving ordinary differential equations I: Nonstiff problems,
    Berlin, New York: Springer-Verlag, ISBN 978-3-540-56670-0

    [2] Hairer, Ernst; Wanner, Gerhard,
    Solving ordinary differential equations II: Stiff and Differential-Algebraic problems,
    Berlin, New York: Springer-Verlag, ISBN 978-3-642-05220-0

  */

  static constexpr int s = 6;
  static constexpr int s_prime = 7;
  
  static constexpr std::array<std::array<double, s>, s>
  a = {0.,            0.,           0.,            0.,         0.,           0.,
       1./5.,         0.,           0.,            0.,         0.,           0.,
       3./40.,        9./40.,       0.,            0.,         0.,           0.,
       44./45.,      -56./15.,      32./9.,        0.,         0.,           0.,
       19372./6561., -25360./2187., 64448./6561., -212./729.,  0.,           0.,
       9017./3168.,  -355./33.,     46732./5247.,  49./176.,  -5103./18656., 0.};
  
  static constexpr std::array<double, s_prime> b = {0., 1./5., 3./10., 4./5., 8./9., 1., 1.};
  static constexpr std::array<double, s> c = {35./384., 0., 500./1113., 125./192., -2187./6784., 11./84.};
  static constexpr std::array<double, s_prime> c_prime = {5179./57600., 0., 7571./16695., 393./640., -92097./339200., 187./2100., 1./40.};

  std::vector<std::vector<double>> k;
  std::vector<double> y_lo;
  std::vector<double> y_hi;
  double err;
  double err_old;
  
  void
  step (double x0,
	std::vector<double>::const_iterator starty0,
	std::vector<double>::const_iterator endy0,
	double dx, ratefun fun) {

    auto ncomp = endy0 - starty0;
    
    static int ii = 0, jj = 0;
    static std::vector<double> ys;    
    k.resize (s_prime);
    for (auto & ik : k)
      ik.assign (ncomp, 0.);
    
    ys.assign (ncomp, 0.);
    y_lo.assign (ncomp, 0.);
    y_hi.assign (ncomp, 0.);
    
    fun (x0, starty0, endy0, k[0].begin ());
    for (ii = 1; ii < s_prime; ++ii) {
      std::copy (starty0, endy0, ys.begin ());
      for (jj = 0; jj < ii; ++jj) {
	auto itk = k[jj].begin ();
	for (auto & refys : ys)
	  refys += *(itk++) * a[ii][jj] * dx;
	//ys += k[jj] * a[ii][jj] * dx;
      }
      fun (x0 + dx * b[ii], ys.begin (), ys.end (), k[ii].begin ());
      //k[ii] = fun (x0 + dx * b[ii], ys);
    }

    std::copy (starty0, endy0, y_lo.begin ());
    for (ii = 0; ii < s; ++ii) {
      auto itk = k[ii].begin ();
      for (auto & refy_lo : y_lo)
	refy_lo += dx * c[ii] * *(itk++);
    }

    std::copy (starty0, endy0, y_hi.begin ());
    for (ii = 0; ii < s_prime; ++ii) {
      auto itk = k[ii].begin ();
      for (auto & refy_hi : y_hi)
	refy_hi += dx * c_prime[ii] * *(itk++);
    }

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
integrate_adaptive (ratefun fun,
		    const std::vector<double> &y0,
		    double xstart, double xend, double tol,
		    std::vector<double> &x, std::vector<double> &y) {

  dorpri5 rk{};
  rk.err_old = tol;
  
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

    rk.step (x.back (), y_back, y.end (), dt, fun);
    

    auto tmpfun  = [] (double x, double y) { return std::abs (x) < std::abs (y); };
    auto normylo = std::abs (*std::max_element (rk.y_lo.begin (), rk.y_lo.end (), tmpfun));

    std::cerr << "t = " << tnew << " dt = " << dt << " dtmin " << dtmin << std::endl;
    if (rk.err < tol * (1. + std::abs (normylo)) || dt <= dtmin) {
      x.push_back (tnew);
      y.insert (y.end (), rk.y_lo.begin (), rk.y_lo.end ());
      y_back = y.end () - ncomp;
      rk.err_old = rk.err;
      //y.push_back (rk.y_lo);
    } else {
      ++rejected;
    }

    //
    // Apply step selection stabilization
    // from [2, formula (2.43c)]
    //
    
    dt = .8 * dt * std::pow (tol / rk.err, 1./5.) * std::pow (rk.err_old / tol, .08);
    
    if (dt < dtmin)
      dt = dtmin;
  }
  std::cerr << rejected << " rejected steps" << std::endl;
};


#endif
