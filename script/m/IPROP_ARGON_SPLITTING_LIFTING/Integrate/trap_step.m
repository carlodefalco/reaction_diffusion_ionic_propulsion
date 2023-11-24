function u1 = trap_step (u0, t, dt, rate, jacrate, mass,nodal_trans)
  b  = mass * u0 + dt * (.001 * rate (u0, t) + .999 * (rate (u0, t) - jacrate (u0, t) * u0)) ;
  A  = mass - .999 * dt * jacrate (u0, t);

  A(nodal_trans,:)=0;
  A(:, nodal_trans)=0;
  A(nodal_trans, nodal_trans)=1;
  b(nodal_trans)=0;
  u1 = full(A) \ b;
endfunction

