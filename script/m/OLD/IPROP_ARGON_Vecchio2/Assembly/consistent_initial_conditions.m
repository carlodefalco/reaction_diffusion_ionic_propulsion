function [u0, udot0_full]= consistent_initial_conditions(y0,t0,x,Vth,mobility,valence, epsilon,q,reactions,index, bc,bc_phi)
  V0t = bc_phi(t0);
  N_species=numfields(index);
  N=numel(x);
%===============================================================================
% DEFINITION OF THE UNKNOWN
%===============================================================================
  n0=[];
  for k=1:N_species
    n0=[ n0, y0(1+N*(k-1):N*k)];
  endfor
%===============================================================================
% ASSEMBLY THE MASS MATRIX FOR EACH SPECIES
%===============================================================================
  Mk_no_bc = bim1a_reaction (x, 1, 1);
  Mk = Mk_no_bc (2:end-1, 2:end-1);

%===============================================================================
% COMPUTE PHI
%===============================================================================
  phi = zeros (N, 1);
  phi([1 end])= V0t;
  P = bim1a_laplacian (x, epsilon, 1);
  b=0;
  for k=1:N_species
    var=valence(k)*n0(:,k);
    b+=var;
  endfor
  b=q*b;


  f_phi_no_bc=b(2:end-1);
  f_phi = Mk_no_bc(2:end-1,2:end-1)*f_phi_no_bc - P(2:end-1,1).*V0t(1)  -  P(2:end-1,end).*V0t(2);
  phi(2:end-1) = P(2:end-1, 2:end-1) \ f_phi;

%===============================================================================
% ASSEMBLY THE GLOBAL MATRICES
%===============================================================================
 % Mass matrix_type
  M=zeros((N_species)*(N),(N_species)*(N));
  for k=1:N_species
    M(1+(N)*(k-1):(N)*k, 1+(N)*(k-1):(N)*k)=Mk_no_bc;
  endfor
  A=zeros(N_species*N, N_species*N );
  for k=1:N_species
    A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, -valence(k)*phi/Vth(k));
  endfor
%===============================================================================
% ASSEMBLY THE RIGHT HAND SIDE
%===============================================================================

  chemistry= compute_change_rates2(n0,reactions,index,x) ;


  F_no_bc= zeros(N*N_species,1);
  for k=1:N_species
    F_no_bc(1+N*(k-1):k*N)  = Mk_no_bc* chemistry (:,k);
  endfor

%===============================================================================
% IMPOSITION OF DIRICHLET BOUNDARY CONDITIONS
%===============================================================================
  F_bc= zeros((N-2)*N_species,1);

  for k=1:N_species
    F_bc(1+(N-2)*(k-1):k*(N-2))= F_no_bc(2+N*(k-1):(k*N)-1) - A(2+N*(k-1):N*k-1, N*(k-1)+1).*bc(k,1) - A(2+N*(k-1):N*k-1,N*k).*bc(k,2);
  endfor

%===============================================================================
% IMPLICIT FUNCTION
%===============================================================================
  u0_bc=[];
  for k=1:N_species
    u0_bc=[u0_bc;n0(2:end-1,k)];
  endfor
  A_bc=zeros(N_species*(N-2),N_species*(N-2));

  M_bc=zeros(N_species*(N-2),N_species*(N-2));
  for k=1:N_species
    A_bc(1+(N-2)*(k-1):(N-2)*k,1+(N-2)*(k-1):(N-2)*k) = A(2+N*(k-1):N*k-1,2+N*(k-1):N*k-1);

    M_bc(1+(N-2)*(k-1):(N-2)*k,1+(N-2)*(k-1):(N-2)*k) = M(2+N*(k-1):N*k-1,2+N*(k-1):N*k-1);
  endfor
  udot0=M_bc\[F_bc-A_bc*u0_bc];
  udot0=[udot0; zeros(N-2,1)];
  BCdot=zeros(N_species+1,2);
  for k=1:N_species
    BCdot(k,1)=F_no_bc(1+(k-1)*N)-Mk_no_bc(1,1)\(A(1+(k-1)*N,1+(k-1)*N)*bc(k,1));
    BCdot(k,2)=F_no_bc(N+(k-1)*N)-Mk_no_bc(end,end)\(A(N+(k-1)*N,N+(k-1)*N)*bc(k,2));
  endfor

  udot0_full=zeros(size(y0));

  for k=1:N_species+1
    udot0_full(1+N*(k-1):k*N)=[BCdot(k, 1); udot0(1+(N-2)*(k-1):k*(N-2));BCdot(k, 2)];
   endfor
   u0=[y0; phi];


 endfunction
