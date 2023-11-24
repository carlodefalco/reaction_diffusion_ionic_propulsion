function u = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass,nodal, boundary_values)

  u(:, 1) = u0;
##  dt = (tspan(2) - tspan(1)) / 100;
dt = ((tspan(end) - tspan(1)) / numel(tspan))/1e1;
  counter=0;
t=tspan(1);
N=size(mass, 1)/7;
  for ispan = 1 : numel (tspan) - 1
    %t  = tspan (ispan)
##    ii = 0;
    while (t <= tspan(ispan + 1))
##      tnew = t + dt;
##      if (tnew >= tspan(ispan + 1))
##        dt   = tspan(ispan + 1) - t;
##	      tnew = tspan(ispan + 1);
##      endif



      unew   = lbwe_step (u0, dt, rate, jacrate, mass,nodal, boundary_values);

      utest1  = lbwe_step (u0, dt/2, rate, jacrate, mass,nodal, boundary_values);

      utest  = lbwe_step (utest1, dt/2, rate, jacrate, mass,nodal, boundary_values);

      %err = norm ((unew - utest) ./ mag, inf)

      err = norm ((unew(end-N+1:end) - utest(end-N+1:end))./mag(end-N+1:end), inf) + norm (log( ( unew(1:end-N)+ones(6*N,1) ) ./ (utest(1:end-N)+ones(6*N,1)) )./ log(mag(1:end-N)), inf);
      err
      t=t+dt;
     u0=unew;
      counter+=1
##      if (err < tol(t))
##        u0 = unew;
##        dt = .6 * (tol(t) / err)^(1/2) * dt
##      else
##
##	      printf ('reject dt_old = %g, dt_new = %g\n', dt, dt/2)
##        dt = dt/2
##      endif

      if(counter==1000)
        counter=0;
        save -binary -z results.gz
      endif
   endwhile
    u(:, ispan+1) = unew;
    error(ispan+1) = err;
  endfor

endfunction
