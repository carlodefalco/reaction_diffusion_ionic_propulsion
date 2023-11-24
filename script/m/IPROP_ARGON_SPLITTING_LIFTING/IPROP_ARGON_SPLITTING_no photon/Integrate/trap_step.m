function u1 = trap_step (u0, t, dt, rate, jacrate, mass)
  b  = mass * u0 + dt * (.001 * rate (u0, t) + .999 * (rate (u0, t) - jacrate (u0, t) * u0)) ;
  A  = mass - .999 * dt * jacrate (u0, t);

  u1 = full(A) \ b;
endfunction

