function udot0= consistent_initial_conditions_ADIM(y0,t0,x,Vth,mobility,valence, epsilon,reactions,index, bc,bc_phi, r_bar, n_bar)
  V0t = bc_phi(t0);
  N_species=numfields(index);
  N=length(x);
%===============================================================================
% DEFINITION OF THE UNKNOWN
%===============================================================================
  n10 = y0(1:N);
  n20 = y0(N+1:2*N);
  n30 = y0(2*N+1:3*N);
  n40 = y0(3*N+1:4*N);
  n50 = y0(4*N+1:5*N);
  n60 = y0(5*N+1:6*N);


  n0=[n10, n20, n30, n40, n50, n60];
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
  b=(valence(1)*n10+valence(2)*n20+valence(3)*n30+valence(4)*n40+valence(5)*n50+valence(6)*n60);
  f_phi_no_bc=b(2:end-1);
  f_phi = Mk_no_bc(2:end-1,2:end-1)*f_phi_no_bc - P(2:end-1,1).*V0t(1)  -  P(2:end-1,end).*V0t(2);
  F_phi= f_phi;
  phi(2:end-1) = P(2:end-1, 2:end-1) \ F_phi;
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

  n0_dim=n0.*n_bar;
  chemistry=[];
  for ii=2:N-1
    chemistry=[chemistry; compute_change_rates(n0_dim(ii,:),reactions,index)] ;
  endfor

  f_1_no_bc   = chemistry (:,1)/r_bar;
  f_2_no_bc   = chemistry (:,2)/r_bar;
  f_3_no_bc   = chemistry (:,3)/r_bar;
  f_4_no_bc   = chemistry (:,4)/r_bar;
  f_5_no_bc   = chemistry (:,5)/r_bar;
  f_6_no_bc   = chemistry (:,6)/r_bar;


%===============================================================================
% IMPOSITION OF DIRICHLET BOUNDARY CONDITIONS
%===============================================================================
  f_1   = f_1_no_bc   - A1(2:end-1,1).*bc(1,1) - A1(2:end-1,end).*bc(1,2);
  f_2   = f_2_no_bc   - A2(2:end-1,1).*bc(2,1) - A2(2:end-1,end).*bc(2,2);
  f_3   = f_3_no_bc   - A3(2:end-1,1).*bc(3,1) - A3(2:end-1,end).*bc(3,2);
  f_4   = f_4_no_bc   - A4(2:end-1,1).*bc(4,1) - A4(2:end-1,end).*bc(4,2);
  f_5   = f_5_no_bc   - A5(2:end-1,1).*bc(5,1) - A5(2:end-1,end).*bc(5,2);
  f_6   = f_6_no_bc
%===============================================================================
% ASSEMBLY THE GLOBAL MATRICES
%===============================================================================
  % Mass matrix_type
  M=zeros((N_species)*(N-2),(N_species)*(N-2));
  Mk=Mk_no_bc(2:end-1, 2:end-1);
  for k=1:N_species
    M(1+(N-2)*(k-1):(N-2)*k, 1+(N-2)*(k-1):(N-2)*k)=Mk;
  endfor

  % Stiffness matrix
  A = zeros((N_species)*(N-2), (N_species)*(N-2));
  A(1:(N-2), 1:(N-2))             = A1(2:end-1,2:end-1);
  A(1+(N-2):2*(N-2), 1+(N-2):2*(N-2))     = A2(2:end-1, 2:end-1);
  A(1+2*(N-2):3*(N-2), 1+2*(N-2):3*(N-2)) = A3(2:end-1, 2:end-1);
  A(1+3*(N-2):4*(N-2), 1+3*(N-2):4*(N-2)) = A4(2:end-1, 2:end-1);
  A(1+4*(N-2):5*(N-2), 1+4*(N-2):5*(N-2)) = A5(2:end-1, 2:end-1);



  F = [Mk*f_1; Mk*f_2; Mk*f_3; Mk*f_4; Mk*f_5; Mk*f_6];

%===============================================================================
% IMPLICIT FUNCTION
%===============================================================================
  u0_bc=[n10(2:end-1); n20(2:end-1); n30(2:end-1); n40(2:end-1); n50(2:end-1); n60(2:end-1)];

  udot0_bc=M\(F-A*u0_bc);


  BCdot1=[f_1_no_bc(1)-Mk_no_bc(1,1)\(A1(1,1)*bc(1,1));
          f_1_no_bc(end)-Mk_no_bc(end,end)\(A1(end,end)*bc(1,2))];
  BCdot2=[f_2_no_bc(1)-Mk_no_bc(1,1)\(A2(1,1)*bc(2,1));
          f_2_no_bc(end)-Mk_no_bc(end,end)\(A2(end,end)*bc(2,2))];
  BCdot3=[f_3_no_bc(1)-Mk_no_bc(1,1)\(A3(1,1)*bc(3,1));
          f_3_no_bc(end)-Mk_no_bc(end,end)\(A3(end,end)*bc(3,2))];
  BCdot4=[f_4_no_bc(1)-Mk_no_bc(1,1)\(A4(1,1)*bc(4,1));
          f_4_no_bc(end)-Mk_no_bc(end,end)\(A4(end,end)*bc(4,2))];
  BCdot5=[f_5_no_bc(1)-Mk_no_bc(1,1)\(A5(1,1)*bc(5,1));
          f_5_no_bc(end)-Mk_no_bc(end,end)\(A5(end,end)*bc(5,2))];
  BCdot6=[f_6_no_bc(1);
          f_6_no_bc(end)];

  udot0=[BCdot1(1); udot0_bc(1:N-2); BCdot1(2);
         BCdot2(1); udot0_bc(1+N-2:(N-2)*2); BCdot2(2);
         BCdot3(1); udot0_bc(1+(N-2)*2:(N-2)*3); BCdot3(2);
         BCdot4(1); udot0_bc(1+(N-2)*3:(N-2)*4); BCdot4(2);
         BCdot5(1); udot0_bc(1+(N-2)*4:(N-2)*5); BCdot5(2);
         BCdot6(1); udot0_bc(1+(N-2)*5:(N-2)*6); BCdot6(2);
         ];


  u0=[n10; n20; n30; n40; n50; n60];
  u0=[u0;phi];
##  u0_bc=[n10(2:end-1); n20(2:end-1); n30(2:end-1); n40(2:end-1); n50(2:end-1); n60(2:end-1)];
##
##  udot0=M\(F-A*u0_bc);
##
##
##  u0=u0_bc


 endfunction
