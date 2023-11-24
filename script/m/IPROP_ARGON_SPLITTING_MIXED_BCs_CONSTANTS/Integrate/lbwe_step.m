function u1 = lbwe_step (u0, dt, rate, jacrate, mass, nodal, boundary_values)

  b  = mass * u0 + dt * (rate (u0) - jacrate (u0) * u0);
  A  = mass - dt * jacrate (u0);

  Abc=A;
  b_bc=b;


%%Robin
  Abc(nodal{2}, nodal{2}-1)= A(nodal{2}, nodal{2}-1)/dt;
  Abc(nodal{2}, nodal{2})= (A(nodal{2}, nodal{2}) -mass(nodal{2}, nodal{2}))/dt + boundary_values{2};%A(nodal{2}, nodal{2})+
  b_bc(nodal{2})=0;
%%Dirichlet
  Abc(nodal{1},:)=[];
  Abc(:, nodal{1})=[];
  b_bc(nodal{1})=[];
  b_bc-=A(setdiff(1:end,nodal{1}),nodal{1})*boundary_values{1};%-A*boundary_values{1};
  warning ('off', 'Octave:nearly-singular-matrix')
  warning('off', 'Octave:singular-matrix')

  u1nobc=Abc\b_bc;
  u1=zeros(size(u0));
  u1(nodal{1})=boundary_values{1};
  keyboard
  u1(setdiff(1:end,nodal{1}))=u1nobc;
 keyboard

endfunction
%% Dirichlet
##  for k=1:numel(nodal{1})
##     Abc(nodal{1}(k),:)=0;
##     Abc(:, nodal{1}(k))=0;
##     Abc(nodal{1}(k), nodal{1}(k))=1;
##     b_bc(nodal{1}(k))=0;
##  endfor
##  u1(nodal{1})=boundary_values{2}(nodal{1});
##  u1 = Abc \ b_bc + boundary_values{1};
