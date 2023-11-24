function [u1,ut_chem, err] = lbwe_step_splitting (u0, t, dt, rate, jacrate, mass, mag)
  N_chem=size(mass{2},1);
  N_trans=size(mass{1},1);
  N=ceil(N_chem/6);
##  N_chem-=2;
  u0_chem=u0(1:N_chem); %qui vanno le componenti del sistema chimico, tutte escluso phi

##  res = @(t, u, udot) mass{2} * udot - rate{2} (u);
##  jac = @(t, u, udot) templatejac (@(t, u, udot) - jacrate{2}(u), @(t, u, udot)mass{2}, t, u, udot);
##  %o = odeset ('InitialStep', dt/100,'Jacobian', jac, 'MaxOrder', 1); %'InitialStep', dt
##  o = odeset ('Jacobian', jac);% 'InitialStep', dt
##  [~, ut_chem] = ode15i (res,[t, t+dt], u0_chem, mass{2} \ rate{2} (u0), o);
##  ut_chem = ut_chem(end, :).';
##
##  if (any(ut_chem) < 0)
##    u1 = u0;
##    err = inf;
##    return
##  endif
ut_chem=u0_chem;
##  ut=[ut_chem(1:N_chem-N); u0(N_chem+1:end)];
  u0(N_chem-N+1:N_chem)=[];
  ut=u0;
  uh  = lbwe_step (ut, dt/2, rate{1}, jacrate{1}, mass{1});
  uh2 = lbwe_step (uh, dt/2, rate{1}, jacrate{1}, mass{1});
  if (any (uh2(1:N_chem)) < 0)
    u1 = u0;
    err = inf;
    return
  endif
  u1  = lbwe_step (ut, dt, rate{1}, jacrate{1}, mass{1});
  err = norm ((uh2(N_chem-N+1:end) - u1(N_chem-N+1:end))./mag(N_chem-N+1:end), inf) + norm (log( ( uh2(1:N_chem-N)+ones(N_chem-N,1) ) ./ (u1(1:N_chem-N)+ones(N_chem-N,1)) )./ log(mag(1:N_chem-N)), inf);
  u1=[u1(1:N_chem-N); ut_chem(N_chem-N+1: N_chem);u1(N_chem-N+1:N_chem)];
endfunction
##  b  = mass{1} * ut + dt * (rate{1} (ut,t) - jacrate{1} (ut,t) * ut) ;
##  A  = mass{1} - dt * jacrate{1} (ut,t);
##  u1 = A \ b;
##
##  b  = mass{1} * ut + dt/2 * (rate{1} (ut,t) - jacrate{1} (ut,t) * ut) ;
##  A  = mass{1} - dt/2 * jacrate{1} (ut,t);
##  uh = A \ b;
##
##  b  = mass{1} * uh + dt/2 * (rate{1} (uh,t) - jacrate{1} (ut,t) * uh) ;
##  A  = mass{1} - dt/2 * jacrate{1} (uh,t);
##  err = norm ((A \ b - u1)./mag, inf)
