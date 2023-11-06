function u1 = lbwe_step (u0, t, dt, rate, jacrate, mass)
  b  = mass * u0 + dt * (rate (u0, t) - jacrate (u0, t) * u0) ;
  A  = mass - dt * jacrate (u0, t);
  u1 = full(A) \ b;
endfunction
