function u1 = lbwe_step (u0, dt, rate, jacrate, mass)

  b  = mass * u0 + dt * (rate (u0) - jacrate (u0) * u0);
  A  = mass - dt * jacrate (u0);
  u1 = A \ b;

  ## for ii = 1:20
  ##   b  = mass * u0 + dt * (rate (u1) - jacrate (u1) * u0) ;
  ##   A  = mass - dt * jacrate (u1);
  ##   unew = A \ b;
  ##   if (norm (unew(1:300)-u1(1:300), inf) < 1e-6);
  ##     u1 = unew;
  ##     break;
  ##   else
  ##     u1 = unew;
  ##   endif
  ## endfor
endfunction

