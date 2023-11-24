function [u, error] = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass, nodal_trans, boundary_values,x,electron_flux)

  u(:, 1) = u0;
  u_chemical(:,1)=u0(1:rows(mass{2}));
  dt = (tspan(2) - tspan(1)) / 100;
##  dt=(tspan(end)-tspan(1))/numel(tspan);
##dt=1e-10;
  counter=0;
  flag=0;
  for ispan = 1 : numel (tspan) - 1
    if ispan==1
      save -binary -z results.gz
    endif
    t  = tspan (ispan)
    ii = 0;

    while (t < tspan(ispan + 1))
      tnew = t + dt;
      if (tnew >= tspan(ispan + 1))
        dt   = tspan(ispan + 1) - t;
	      tnew = tspan(ispan + 1);
      endif


      [unew,ut_chem, err]   = lbwe_step_splitting (u0, t, dt, rate, jacrate, mass, mag, nodal_trans, boundary_values,x,electron_flux);
      counter+=1
      err
      if (err < tol(t)  )%&&t/dt<1e3
        t=tnew
        u0 = unew;
        dt = .6 * (tol(t) / err)^(1/2) * dt
##      elseif (t/dt>1e3)
##        t=tnew
##        u0=unew;
##        dt=t/100;
##        flag=1;
##      else
##        if flag==1
##          t=tnew
##          u0=unew;
##          dt=t/100;
        else
        printf ('reject dt_old = %g, dt_new = %g\n', dt, dt/2)
        dt = dt/2
####         endif
        endif

    if(counter==10)
      counter=0;
      save -binary -z results.gz
    endif
    endwhile
    u_chemical(:, ispan+1)=ut_chem;
    u(:, ispan+1) = unew;
    error(ispan+1) = err;

endfor
save -binary -z results.gz
endfunction


