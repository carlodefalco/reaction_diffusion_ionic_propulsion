function [u,error] = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass,r,idx,x,N_species,database,nodal_trans,boundary_values)

  u(:, 1) = u0;
  u_chemical(:,1)=u0(1:rows(mass{2}));
  Te(:,1)= compute_temperature( u0,x,idx,N_species, database);
##  dt = (tspan(2) - tspan(1)) / 100;
dt=(tspan(end)-tspan(1))/numel(tspan);
  counter=0;
  flag=0;
  for ispan = 1 : numel (tspan) - 1
    if ispan==1
      save -binary -z results.gz
    endif
    t  = tspan (ispan);
    ii = 0;
##    while (t < tspan(ispan + 1))
##      tnew = t + dt;
##      if (tnew >= tspan(ispan + 1))
##        dt   = tspan(ispan + 1) - t;
##	      tnew = tspan(ispan + 1);
##      endif
      r_new= compute_coeff(u0, x,r, idx,N_species, database);


      [unew,ut_chem, err]   = lbwe_step_splitting (u0, t, dt, rate, jacrate, mass, mag,database,r_new,nodal_trans,boundary_values);
      counter+=1
      err
      u0=unew;
      t+=dt
      dt


##      if (err < tol(t) && t/dt<1e3)
##        t=tnew
##        u0 = unew;
##        dt = .6 * (tol(t) / err)^(1/2) * dt
##      elseif (t/dt>1e3)
##        t=tnew
##        u0=unew;
##        dt=t/100
##        flag=1;
##      else
##        if flag==1
##          t=tnew
##          u0=unew;
##          dt=t/100
##        else
##          printf ('reject dt_old = %g, dt_new = %g\n', dt, dt/2)
##          dt = dt/2
##        endif
##      endif
      %if(counter==50)
      %counter=0;
      save -binary -z results.gz
    %endif
##  endwhile
  u_chemical(:, ispan+1)=ut_chem;
  Te(:, ispan+1)= compute_temperature( unew,x,idx,N_species, database);
  u(:, ispan+1) = unew;
  error(ispan+1) = err;
endfor

endfunction
