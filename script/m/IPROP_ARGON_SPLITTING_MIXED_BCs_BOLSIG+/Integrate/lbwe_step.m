function u1 = lbwe_step (u0, dt, rate, jacrate, mass,nodal_trans,boundary_values)
  b  = mass * u0 + dt * (rate (u0) - jacrate (u0) * u0);
  A  = mass - dt * jacrate (u0);
  Abc=A;

  b_bc=b-A*boundary_values;

  for k=1:numel(nodal_trans)
     Abc(nodal_trans(k),:)=0;
     Abc(:, nodal_trans(k))=0;
     Abc(nodal_trans(k), nodal_trans(k))=1;
     b_bc(nodal_trans(k))=0;
  endfor

  warning ('off', 'Octave:nearly-singular-matrix')
  warning('off', 'Octave:singular-matrix')
  u1 = Abc \ b_bc + boundary_values;

  ## for ii = 1:20
  ##   b  = mass * u0 + dt * (rate (u1) - jacrate (u1) * u0) ;
  ##   A  = mass - dt * jacrate (u1);
  ##   unew = A \ b;
  ##   if (norm (unew(1:300)-u1(1:300), inf) < 1e-6);
  ##     u1 = unew;
  ##     break;
  ##   else
  ##     u1 = unew;
  ##   endif
  ## endfor
endfunction

