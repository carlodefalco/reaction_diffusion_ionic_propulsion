function eqs = compute_drift_diffusion_reaction_system (t, y, ydot, reactions, index, x, N, q, epsilon, Vth, mobility, valence, bc, bc_phi)
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
  N_state=N_species+1;
%===============================================================================
% DEFINITION OF THE UNKNOWN
%===============================================================================

% rearrange the unknown in a matrix without the boundary value
  n=[];
  ndot=[];
  for k=1:N_species
    n    = [n,y(2+(k-1)*N: k*N-1)];  %remove left and right boundary from y
    ndot = [ndot,ydot(2+(k-1)*N: k*N-1)];  %remove left and right boundary from ydot
  endfor
%===============================================================================
% ASSEMBLY THE MASS MATRIX FOR EACH SPECIES
%===============================================================================
  Mk_no_bc = bim1a_reaction (x, 1, 1);


%===============================================================================
% ASSEMBLY THE GLOBAL MATRICES
%===============================================================================
  phi=[V0t(1); n(:,end); V0t(2)];
  A=zeros(N_species*N, N_species*N );
  for k=1:N_species
    A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, -valence(k)*phi/Vth(k));
  endfor
  P = bim1a_laplacian (x, epsilon, 1);


%===============================================================================
% ASSEMBLY THE RIGHT HAND SIDE
%===============================================================================
  b_no_bc=0;
  for k=1:N_species
    var=valence(k)*n(:,k);
    b_no_bc+=var;
  endfor
  b_no_bc=q*b_no_bc;


  chemistry=compute_change_rates2(n,reactions,index,x(2:end-1)) ;

  F_no_bc= zeros((N-2)*N_species,1);
  for k=1:N_species
    F_no_bc(1+(N-2)*(k-1):k*(N-2))  = Mk_no_bc(2:end-1, 2:end-1)* chemistry (:,k);
  endfor

  f_phi_no_bc = b_no_bc;
keyboard
%===============================================================================
%IMPOSITION OF THE BOUNDARY CONDITIONS
%===============================================================================
 Mk=Mk_no_bc(2:end-1, 2:end-1);

  f_phi = Mk*f_phi_no_bc - P(2:end-1,1).*V0t(1)  -  P(2:end-1,end).*V0t(2);
  F_bc= zeros((N-2)*N_species,1);

  for k=1:N_species
    F_bc(1+(N-2)*(k-1):k*(N-2))= F_no_bc(1+(N-2)*(k-1):(k*(N-2))) - A(2+N*(k-1):N*k-1, N*(k-1)+1).*bc(k,1) - A(2+N*(k-1):N*k-1,N*k).*bc(k,2);
  endfor


%===============================================================================
% IMPLICIT FUNCTION
%===============================================================================
##  u=[n1; n2; n3; n4; n5; n6; phi];
##  udot=[n1dot; n2dot; n3dot; n4dot; n5dot; n6dot; zeros(N-2,1)];
##  eqs=M*udot+A*u-F;
 M=zeros((N_species)*(N-2),(N_species)*(N-2));
  for k=1:N_species
    M(1+(N-2)*(k-1):(N-2)*k, 1+(N-2)*(k-1):(N-2)*k)=Mk;
  endfor

  A_bc=zeros(N_species*(N-2),N_species*(N-2));


  for k=1:N_species
    A_bc(1+(N-2)*(k-1):(N-2)*k,1+(N-2)*(k-1):(N-2)*k) = A(2+N*(k-1):N*k-1,2+N*(k-1):N*k-1);
  endfor

  y_bc=[];
  ydot_bc=[];
  for k=1:N_species
    y_bc    = [y_bc; n(:, k)];
    ydot_bc = [ydot_bc; ndot(:,k)];
  endfor

keyboard
 eqs= [M*ydot_bc+A_bc*y_bc-F_bc;
        P(2:end-1, 2:end-1)*y(N_species*N+2:end-1)-f_phi];
keyboard
 endfunction

