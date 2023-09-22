function [u1,err] = lbwe_step_splitting (u0, t, dt, rate, jacrate, mass, mag)
  N_chem=size(mass{2},1);
  u0_chem=u0(1:N_chem); %qui vanno le componenti del sistema chimico, tutte escluso phi

  res = @(t, u, udot) mass{2} * udot - rate{2} (u,t);
  jac = @(t, u, udot) templatejac (@(t, u, udot) - jacrate{2}(u,t), @(t, u, udot)mass{2}, t, u, udot);
  o = odeset ('Jacobian', jac, 'InitialStep', dt, 'MaxOrder', 1);
  [~, ut_chem] = ode15i (res, [t, t+dt], u0_chem, mass{2} \ rate{2} (u0,t), o);
  ut_chem = ut_chem(end, :).';
  ut=[ut_chem; u0(N_chem+1:end)];

  b  = mass{1} * ut + dt * (rate{1} (ut,t) - jacrate{1} (ut,t) * ut) ;
  A  = mass{1} - dt * jacrate{1} (ut,t);
  u1 = A \ b;

  b  = mass{1} * ut + dt/2 * (rate{1} (ut,t) - jacrate{1} (ut,t) * ut) ;
  A  = mass{1} - dt/2 * jacrate{1} (ut,t);
  uh = A \ b;

  b  = mass{1} * uh + dt/2 * (rate{1} (uh,t) - jacrate{1} (ut,t) * uh) ;
  A  = mass{1} - dt/2 * jacrate{1} (uh,t);
  err = norm ((A \ b - u1)./mag, inf);

endfunction
