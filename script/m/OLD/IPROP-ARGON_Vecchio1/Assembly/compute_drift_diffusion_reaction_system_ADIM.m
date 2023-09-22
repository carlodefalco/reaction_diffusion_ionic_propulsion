function eqs = compute_drift_diffusion_reaction_system_ADIM (t, y, ydot, reactions, index, x, N, epsilon, Vth, mobility, valence, bc, bc_phi, r_bar,n_bar)
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
%          x             : (double)              array of grid points
%          N             : (int)                 number of grid points
%          V0            : (function handle)     applied tension
%          q             : (double)              fundamental charge constant
%          epsilon       : (double)              dielectric constant
%          Vth           : (double)
%          mobility      : (double)              mobility of the species
%          valence       : (double)              valence number of the species
%          bc            : (double)              matrix of boundary conditions
%    OUTPUT:
%          eqs          : (function handle)     function to be passed to ode15i
% Usage:
%    eqs =@(t, y, ydot) compute_drift_diffusion_reaction_system (t, y, ydot,
%      reactions,index, x, N, V0, q, epsilon, Vth, mobility, valence, bc)
%   [t, sol]=ode15i(eqs, ...)

  V0t=bc_phi(t);
  N_species=numfields(index);

%===============================================================================
% DEFINITION OF THE UNKNOWN
%===============================================================================
  n1 = y(1:(N-2));
  n2 = y((N-2)+1:2*(N-2));
  n3 = y(2*(N-2)+1:3*(N-2));
  n4 = y(3*(N-2)+1:4*(N-2));
  n5 = y(4*(N-2)+1:5*(N-2));
  n6 = y(5*(N-2)+1:6*(N-2));

  n1dot=ydot(1:(N-2));
  n2dot=ydot((N-2)+1:2*(N-2));
  n3dot=ydot(2*(N-2)+1:3*(N-2));
  n4dot=ydot(3*(N-2)+1:4*(N-2));
  n5dot=ydot(4*(N-2)+1:5*(N-2));
  n6dot=ydot(5*(N-2)+1:6*(N-2));

  n=[n1, n2, n3, n4, n5, n6];
%===============================================================================
% ASSEMBLY THE MASS MATRIX FOR EACH SPECIES
%===============================================================================
  Mk_no_bc = bim1a_reaction (x, 1, 1);
%===============================================================================
% COMPUTE PHI
%===============================================================================
  phi = zeros (N, 1);
  phi([1 end])= V0t;
  P = bim1a_laplacian (x, epsilon, 1);
  b=(valence(1)*n1+valence(2)*n2+valence(3)*n3+valence(4)*n4+valence(5)*n5+valence(6)*n6);
  f_phi_no_bc=b;
  f_phi = Mk_no_bc(2:end-1,2:end-1)*f_phi_no_bc - P(2:end-1,1).*V0t(1)  -  P(2:end-1,end).*V0t(2);
  F_phi= f_phi;
  phi(2:end-1) = P(2:end-1, 2:end-1) \ F_phi;
%===============================================================================
%===============================================================================
% ASSEMBLY THE STIFFNESS MATRIX FOR EACH SPECIES
%===============================================================================

  A1 = bim1a_advection_diffusion(x, mobility(1)*Vth(1), 1, 1, -valence(1)*phi/Vth(1));
  %dn1(2:end-1) = A1(2:end-1, :) * n1;

  A2 = bim1a_advection_diffusion(x, mobility(2)*Vth(2), 1, 1, 0);
  %dn2(2:end-1) = A2(2:end-1, :) * n2;

  A3 = bim1a_advection_diffusion(x, mobility(3)*Vth(2), 1, 1, -valence(3)*phi/Vth(2));
  %dn3(2:end-1) = A3(2:end-1, :) * n3;

  A4 = bim1a_advection_diffusion(x, mobility(4)*Vth(2), 1, 1, 0);
  %dn4(2:end-1) = A4(2:end-1, :) * n4;

  A5 = bim1a_advection_diffusion(x, mobility(5)*Vth(2), 1, 1, -valence(5)*phi/Vth(2));
  %dn5(2:end-1) = A5(2:end-1, :) * n5;





%===============================================================================
% ASSEMBLY THE RIGHT HAND SIDE
%===============================================================================


  n_dim=n.*n_bar;
  chemistry=[];
  for ii=1:N-2
    chemistry=[chemistry; compute_change_rates(n_dim(ii,:),reactions,index)] ;
  endfor

  f_1_no_bc   = chemistry (:,1)/r_bar;
  f_2_no_bc   = chemistry (:,2)/r_bar;
  f_3_no_bc   = chemistry (:,3)/r_bar;
  f_4_no_bc   = chemistry (:,4)/r_bar;
  f_5_no_bc   = chemistry (:,5)/r_bar;
  f_6_no_bc   = chemistry (:,6)/r_bar;

%===============================================================================
% ASSEMBLY THE GLOBAL MATRICES
%===============================================================================
  Mk=Mk_no_bc(2:end-1, 2:end-1);
%===============================================================================
%IMPOSITION OF THE BOUNDARY CONDITIONS
%===============================================================================

  f_1   = Mk*f_1_no_bc - A1(2:end-1,1).*bc(1,1) - A1(2:end-1,end).*bc(1,2);
  f_2   = Mk*f_2_no_bc - A2(2:end-1,1).*bc(2,1) - A2(2:end-1,end).*bc(2,2);
  f_3   = Mk*f_3_no_bc - A3(2:end-1,1).*bc(3,1) - A3(2:end-1,end).*bc(3,2);
  f_4   = Mk*f_4_no_bc - A4(2:end-1,1).*bc(4,1) - A4(2:end-1,end).*bc(4,2);
  f_5   = Mk*f_5_no_bc - A5(2:end-1,1).*bc(5,1) - A5(2:end-1,end).*bc(5,2);
  f_6   = Mk*f_6_no_bc;


%===============================================================================
% IMPLICIT FUNCTION
%===============================================================================
  eqs=[Mk*n1dot+A1(2:end-1,2:end-1)*n1-f_1;
       Mk*n2dot+A2(2:end-1,2:end-1)*n2-f_2;
       Mk*n3dot+A3(2:end-1,2:end-1)*n3-f_3;
       Mk*n4dot+A4(2:end-1,2:end-1)*n4-f_4;
       Mk*n5dot+A5(2:end-1,2:end-1)*n5-f_5;
       Mk*n6dot-f_6];



 endfunction
