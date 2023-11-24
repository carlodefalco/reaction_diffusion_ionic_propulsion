close all;
clear all;


function u1 = lbwe_step (u0, t, dt, rate, jacrate, mass)
  b  = mass * u0 + dt * (rate (u0, t) - jacrate (u0, t) * u0) ;
  A  = mass - dt * jacrate (u0, t);
  u1 = A \ b;
endfunction

function u = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass)

  u(:, 1) = u0;
  dt = (tspan(2) - tspan(1)) / 10;

  h = waitbar (0);
  for ispan = 1 : numel (tspan) - 1
    t  = tspan (ispan);
    while (t < tspan(ispan + 1))
      if (t + dt > tspan(ispan + 1) - 100*eps)
        dt = tspan(ispan + 1)  - t;
      endif
      t = t + dt;

      unew   = lbwe_step (u0, t, dt, rate, jacrate, mass);

      utest  = lbwe_step (u0, t-dt/2, dt/2, rate, jacrate, mass);
      utest  = lbwe_step (utest, t, dt/2, rate, jacrate, mass);

      err = norm (unew - utest, inf) / mag;

      if (err < tol)
        u0 = unew;
        dt = .6 * sqrt (tol / err) * dt
        h = waitbar ((t-tspan(1))/(tspan(end)-tspan(1)), h);
      else
        t = t-dt;
        dt = min (.6 * sqrt (tol / err), .5) * dt;
      endif
      ispan
    endwhile
    u(:, ispan+1) = unew;
    endfor

endfunction


u0 = [1; 1];
tspan = (linspace (0, 10, 100));
tol = 1e-7;
mag = norm (u0, inf);
rate = @(u, t) [-3*u(1); -2*u(2)];
jacrate = @(u, t) [-3, 0; 0, -2];
mass = [1, 0; 0, 1];

u = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass);
semilogy (tspan, u, tspan, [exp(-3*tspan); exp(-2*tspan)], 'x')


