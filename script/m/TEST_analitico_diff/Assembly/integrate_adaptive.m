function u = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass)

  u(:, 1) = u0;
  dt = (tspan(2) - tspan(1)) / 100;

  %%h = waitbar (0);
  for ispan = 1 : numel (tspan) - 1
    t  = tspan (ispan);
    ii = 0;
    while (t < tspan(ispan + 1))
      if (t + dt > tspan(ispan + 1) - 100*eps)
	tnew = tspan(ispan + 1);
        dt = tnew - t;
	t  = tnew;
      endif
      t = t + dt

      unew   = trap_step (u0, t, dt, rate, jacrate, mass);

      utest1  = trap_step (u0, t-dt/2, dt/2, rate, jacrate, mass);

      utest  = trap_step (utest1, t, dt/2, rate, jacrate, mass);

      err = norm ((unew - utest) ./ mag, inf)

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
##    N = 300;
##    L = 1.5e-6;
##    x = linspace (0, L, N).';
##    plot (x(2:end-1),  (u(1:2*N, [ispan ispan+1]), 1, 2),
##	  x, diff (u(1: N-2, [ispan ispan+1]), 1, 2),
##	  x, mag(1: N-2) * tol(t) .* [-1, 1])
##    drawnow
  endfor
  %%close (h)
endfunction
