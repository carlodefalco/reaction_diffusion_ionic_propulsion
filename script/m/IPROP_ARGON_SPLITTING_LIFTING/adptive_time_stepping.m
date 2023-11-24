close all;
clear all;


function [dfdy, dfdydot] = templatejac (dy, dydot, t, y, ydot)
    dfdy = dy (t, y, ydot) ;
    if (nargout>1)
      dfdydot = dydot (t, y, ydot) ;
    endif
endfunction

function [u1, err] = lbwe_step_splitting (u0, t, dt, rate, jacrate, mass, mag)

  res = @(t, u, udot) mass * udot - rate{2} (u);
  jac = @(t, u, udot) templatejac (@(t, u, udot) - jacrate{2} (u), @(t, u, udot) mass, t, u, udot);
  o = odeset ('Jacobian', jac, 'InitialStep', dt, 'MaxOrder', 1);
  [~, ut] = ode15i (res, [t, t+dt], u0, mass \ rate{2} (u0), o);
  rows (ut)
  ut = ut(end, :).';

  b  = mass * ut + dt * (rate{1} (ut) - jacrate{1} (ut) * ut) ;
  A  = mass - dt * jacrate{1} (ut);
  u1 = A \ b;

  b  = mass * ut + dt/2 * (rate{1} (ut) - jacrate{1} (ut) * ut) ;
  A  = mass - dt/2 * jacrate{1} (ut);
  uh = A \ b;

  b  = mass * uh + dt/2 * (rate{1} (uh) - jacrate{1} (ut) * uh) ;
  A  = mass - dt/2 * jacrate{1} (uh);
  err = norm ((A \ b - u1)./mag, inf);

endfunction

%%mass (u1 - u0) / dt = 1/2 (rate (u0, t) )

function u = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass)

  u(:, 1) = u0;
  dt = (tspan(2) - tspan(1)) / 100;

  %%h = waitbar (0);
  for ispan = 1 : numel (tspan) - 1
    t  = tspan (ispan);
    ii = 0;
    while (t < tspan(ispan + 1))
      tnew = t + dt;
      if (tnew >= tspan(ispan + 1))
        dt   = tspan(ispan + 1) - t;
	      tnew = tspan(ispan + 1);
      endif
      t = tnew

      [unew, err]   = lbwe_step_splitting (u0, t, dt, rate, jacrate, mass, mag);

      ## utest  = lbwe_step_splitting (u0, t-dt/2, dt/2, rate, jacrate, mass);
      ## utest  = lbwe_step_splitting (utest, t, dt/2, rate, jacrate, mass);
      ## err = norm ((unew - utest) ./ mag, inf)

      if (err < tol(t))
        u0 = unew;
        dt = .6 * (tol(t) / err)^(1/2) * dt
      else
        t = t-dt;
	printf ('reject dt_old = %g, dt_new = %g\n', dt, dt/2)
        dt = dt/2
      endif

    endwhile
    u(:, ispan+1) = unew;
  endfor
  %%close (h)
endfunction


u0 = [1; 1];
rate = {@(u) [-3*u(1); -2*u(2)], @(u) [-1*u(1); -2*u(2)]};
jacrate = {@(u) [-3, 0; 0, -2], @(u) [-1, 0; 0, -2]};
mass = eye(2);
tol  = @(t) exp (-4 * t);
mag = .01;
tspan = linspace (0, 10, 10);
u = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass);
semilogy (tspan, u, tspan, exp (-4*tspan), '-.')
'x')


