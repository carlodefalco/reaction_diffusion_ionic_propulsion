function eqs = compute_drift_diffusion_reaction_system_2 (t, y, ydot, reactions, index, x, N, V0, q, epsilon, Vth, mobility, valence)
%===============================================================================
% Assembly of the system for simulating the plasma discharge with finite elements
% (homogeneous Dirichlet boundary conditions) to used with ode15i.
% The system is build such that R(t, u, u')=0;
% Mu'+Au=F---> R(t,u, u')=Mu'+Au-F; u=[n1,n2...n_N, phi]'
%===============================================================================
%
%    INPUT:
%          t             : (double)              time variable
%          y             : (double)              unknown
%          ydot          : (double)              time derivative of the unknown
%          reactions     : (struct)              reactions considered
%          idex          : (struct)              index of each species
%          x             : (int)                 array of grid points
%          N             : (int)                 number of grid points
%          V0            : (function handle)     applied tension
%          q             : (int)                 fundamental charge constant
%          epsilon       : (int)                 dielectric constant
%          Vth_e         : (int)
%          Vth_i         : (int)
%          mobility      : (int)                 mobility of the species
%          valence       : (double)              valence number of the species
%    OUTPUT:
%          eqs          : (function handle)     function to be passed to ode15i
% Usage:
%    eqs =@(t, y, ydot) compute_drift_diffusion_reaction_system (t, y, ydot,
%      reactions,index, x, N, V0, q, epsilon, Vths, mobility, valence)
%   [t, sol]=ode15i(eqs, ...)

  V0t = V0(t);
  N_species=numfields(index);

%===============================================================================
% DEFINITION OF THE UNKNOWN
%===============================================================================
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
  phidot= ydot(6*N+1:7*N);

  n=[n1, n2, n3, n4, n5, n6];

%===============================================================================
% ASSEMBLY THE STIFFNESS MATRIX FOR EACH SPECIES
%===============================================================================

  A1 =- bim1a_advection_diffusion(x, mobility(1)*Vth(1), 1, 1, -valence(1)*phi/Vth(1));
  %dn1(2:end-1) = A1(2:end-1, :) * n1;

  A2 = -bim1a_advection_diffusion(x, mobility(2)*Vth(2), 1, 1, -valence(2)*phi/Vth(2));
  %dn2(2:end-1) = A2(2:end-1, :) * n2;

  A3 = -bim1a_advection_diffusion(x, mobility(3)*Vth(2), 1, 1, -valence(3)*phi/Vth(2));
  %dn3(2:end-1) = A3(2:end-1, :) * n3;

  A4 = -bim1a_advection_diffusion(x, mobility(4)*Vth(2), 1, 1, -valence(4)*phi/Vth(2));
  %dn4(2:end-1) = A4(2:end-1, :) * n4;

  A5 = -bim1a_advection_diffusion(x, mobility(5)*Vth(2), 1, 1, -valence(5)*phi/Vth(2));
  %dn5(2:end-1) = A5(2:end-1, :) * n5;

  A6 = -bim1a_advection_diffusion(x, mobility(6)*Vth(2), 1, 1, -valence(6)*phi/Vth(2));
  %dn6(2:end-1) = A6(2:end-1, :) * n6;

  P = bim1a_laplacian (x, epsilon, 1);
  %phi([1 N]) = V0t;


%===============================================================================
% ASSEMBLY THE MASS MATRIX FOR EACH SPECIES
%===============================================================================
  Mk = bim1a_reaction (x, 1, 1);

%===============================================================================
% IMPOSITION OF DIRICHLET BOUNDARY CONDITIONS
%===============================================================================
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

%===============================================================================
% ASSEMBLY THE GLOBAL MATRICES
%===============================================================================
  % Mass matrix_type
  M=zeros((N_species+1)*N,(N_species+1)*N);
  for k=1:N_species
    M(1+N*(k-1):N*k, 1+N*(k-1):N*k)=Mk;
  endfor

  % Stiffness matrix
  A = zeros(N_species*N+N, N_species*N+N);
  A(1:N, 1:N) = A1;
  A(1+N:2*N, 1+N:2*N)= A2;
  A(1+2*N:3*N, 1+2*N:3*N) = A3;
  A(1+3*N:4*N, 1+3*N:4*N) = A4;
  A(1+4*N:5*N, 1+4*N:5*N) = A5;
  A(1+5*N:6*N, 1+5*N:6*N) = A6;
  A(1+6*N:7*N, 1+6*N:7*N) = P;

%===============================================================================
% ASSEMBLY THE RIGHT HAND SIDE
%===============================================================================
  b=q*(valence(1)*n1+valence(2)*n2+valence(3)*n3+valence(4)*n4+valence(5)*n5+valence(6)*n6);
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

  f_1([1 end]) = [0 0];
  f_2([1 end]) = [0 0];
  f_3([1 end]) = [0 0];
  f_4([1 end]) = [0 0];
  f_5([1 end]) = [0 0];
  f_6([1 end]) = [0 0];
  f_phi        = b;
  f_phi([1 end]) = V0t;


  F = [Mk*f_1; Mk*f_2; Mk*f_3; Mk*f_4; Mk*f_5; Mk*f_6; Mk*f_phi];

  %phi(2:end-1) = A00(2:end-1, 2:end-1) \ (b(2:end-1)-A00(2:end-1, :) * phi);
%===============================================================================
% IMPLICIT FUNCTION
%===============================================================================
  u=[n1; n2; n3; n4; n5; n6; phi];
  udot=[n1dot; n2dot; n3dot; n4dot; n5dot; n6dot; phidot];
  eqs=M*udot+A*u-F;

 endfunction
