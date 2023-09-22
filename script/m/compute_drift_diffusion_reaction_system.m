%% Solve :
%%   n_k' + div (mun Vth (grad n_k - z_k grad phi/Vth n_k)) = ndot
%%   0  = -div (epsilon grad (phi)) - q sum_k { z_k  n_k}
%%   n_k(0) = n_k(L)= 0
%%   phi(0) = 0 phi(L) = V0(t)
%%   this function builds the system of ODE in implicit form

function eqs = compute_drift_diffusion_reaction_system (t, y, ydot, reactions, index, x, N, V0, q, epsilon, Vth, mobility, valence)

  V0t = V0(t);
  N_species=numfields(index);
  n1 = y(1:N);
  n2 = y(N+1:2*N);
  n3 = y(2*N+1:3*N);
  n4 = y(3*N+1:4*N);
  n5 = y(4*N+1:5*N);
  n6 = y(5*N+1:6*N);
  phi= y(6*N+1:7*N);

  n1dot=ydot(1:N);
  n2dot=ydot(N+1:2*N);
  n3dot=ydot(2*N+1:3*N);
  n4dot=ydot(3*N+1:4*N);
  n5dot=ydot(4*N+1:5*N);
  n6dot=ydot(5*N+1:6*N);

  n=[n1, n2, n3, n4, n5, n6];

  Mk = bim1a_reaction (x, 1, 1);
  sum_k=[];
  for k=1:N_species
    sum_k= [sum_k, n(:, k) .* valence(k)];
  endfor
  sum_k=sum(sum_k,2);
  b=q*Mk*sum_k;



  A1 = bim1a_advection_diffusion(x, mobility(1)*Vth, 1, 1, -valence(1)*phi/Vth);
  %dn1(2:end-1) = A1(2:end-1, :) * n1;

  A2 = bim1a_advection_diffusion(x, mobility(2)*Vth, 1, 1, -valence(2)*phi/Vth);
  %dn2(2:end-1) = A2(2:end-1, :) * n2;

  A3 = bim1a_advection_diffusion(x, mobility(3)*Vth, 1, 1, -valence(3)*phi/Vth);
  %dn3(2:end-1) = A3(2:end-1, :) * n3;

  A4 = bim1a_advection_diffusion(x, mobility(4)*Vth, 1, 1, -valence(4)*phi/Vth);
  %dn4(2:end-1) = A4(2:end-1, :) * n4;

  A5 = bim1a_advection_diffusion(x, mobility(5)*Vth, 1, 1, -valence(5)*phi/Vth);
  %dn5(2:end-1) = A5(2:end-1, :) * n5;

  A6 = bim1a_advection_diffusion(x, mobility(6)*Vth, 1, 1, -valence(6)*phi/Vth);
  %dn6(2:end-1) = A6(2:end-1, :) * n6;
  P = bim1a_laplacian (x, epsilon, 1);
  %phi([1 N]) = V0t;

  %imposition of the boundary condtions
  A1([1 end], :) = 0;
  A2([1 end], :) = 0;
  A3([1 end], :) = 0;
  A4([1 end], :) = 0;
  A5([1 end], :) = 0;
  A6([1 end], :) = 0;
  P ([1 end], :) = 0;

  A1(1,1)=1;
  A2(1,1)=1;
  A3(1,1)=1;
  A4(1,1)=1;
  A5(1,1)=1;
  A6(1,1)=1;
  P (1,1)=1;

  A1(end,end)=1;
  A2(end,end)=1;
  A3(end,end)=1;
  A4(end,end)=1;
  A5(end,end)=1;
  A6(end,end)=1;
  P (end,end)=1;
  %assemble the global matrices
  M=zeros((N_species+1)*N,(N_species+1)*N);
  for k=1:N_species
    M(1+N*(k-1):N*k, 1+N*(k-1):N*k)=Mk;
  endfor
  A = zeros(N_species*N+N, N_species*N+N);
  A(1:N, 1:N) = A1;
  A(1+N:2*N, 1+N:2*N)= A2;
  A(1+2*N:3*N, 1+2*N:3*N) = A3;
  A(1+3*N:4*N, 1+2*N:3*N) = A4;
  A(1+4*N:5*N, 1+2*N:3*N) = A5;
  A(1+5*N:6*N, 1+2*N:3*N) = A6;
  A(1+6*N:7*N, 1+2*N:3*N) = P;
  %right hand side
  chemistry=[];
  for ii=1:N
    chemistry=[chemistry; compute_change_rates(n(ii,:),reactions,index)] ;
  endfor
  f_1 = chemistry(:,1);
  f_2 = chemistry(:,2);
  f_3 = chemistry(:,3);
  f_4 = chemistry(:,4);
  f_5 = chemistry(:,5);
  f_6 = chemistry(:,6);

  f_1([1 end])= [0 0];
  f_2([1 end])= [0 0];
  f_3([1 end])= [0 0];
  f_4([1 end])= [0 0];
  f_5([1 end])= [0 0];
  f_6([1 end])= [0 0];
  f_phi=b;
  f_phi([1 end])= V0t;
  u=[n1; n2; n3; n4; n5; n6; phi];
  udot=[n1dot; n2dot; n3dot; n4dot; n5dot; n6dot; zeros(N,1)];

  F=M*[f_1; f_2; f_3; f_4; f_5; f_6; f_phi];

  %phi(2:end-1) = A00(2:end-1, 2:end-1) \ (b(2:end-1)-A00(2:end-1, :) * phi);

  eqs=M*udot+A*u-F;
 endfunction
