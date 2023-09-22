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

##      utest  = lbwe_step_splitting (u0, t-dt/2, dt/2, rate, jacrate, mass);
##      utest  = lbwe_step_splitting (utest, t, dt/2, rate, jacrate, mass);
##      err = norm ((unew - utest) ./ mag, inf)

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
