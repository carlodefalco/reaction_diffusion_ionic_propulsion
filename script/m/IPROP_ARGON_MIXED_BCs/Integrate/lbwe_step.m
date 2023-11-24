function u1 = lbwe_step (u0, dt, rate, jacrate, mass, nodal, boundary_values)
  b  = mass * u0 + dt * (rate (u0) - jacrate (u0) * u0);
  A  = mass - dt * jacrate (u0);
  Abc=A;

  b_bc=b-A*boundary_values;

  for k=1:numel(nodal)
     Abc(nodal(k),:)=0;
     Abc(:, nodal(k))=0;
     Abc(nodal(k), nodal(k))=1;
     b_bc(nodal(k))=0;
  endfor

  warning ('off', 'Octave:nearly-singular-matrix')
  warning('off', 'Octave:singular-matrix')
  u1 = Abc \ b_bc + boundary_values;
endfunction
