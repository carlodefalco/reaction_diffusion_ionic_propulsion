
%% Solve :
%%   n_k' = div (mun Vth (grad n_k - grad phi/Vth n_k))
function [J, DJac]= compute_jacobian_system_2 (t, y, ydot, reactions, index, x, N, V0, q, epsilon, Vth, mobility, valence)
%===============================================================================
% Assembly of the Jacobian of the system for simulating the plasma discharge
% with finite elements(homogeneous Dirichlet boundary conditions) to used with
% the option of ode15i.
% The system is build such that R(t, u, u')=0;
% Mu'+Au=F---> R(t,u, u')=Mu'+Au-F; u=[n1,n2...n_N, phi]'
% The Jacobian is J=dR_du while DJac=dR_du';
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
%          [J,DJac]      : (function handle)     function to be passed to odeset
% Usage:
%   options=odeset('Jacobian', @(t, y, ydot)  compute_jacobian_system (t, y,
%   ydot,reactions,index, x, N, V0, q, epsilon, Vth_e,Vth_ions, mobility, valence))
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

  n=[n1, n2, n3, n4, n5, n6];

%===============================================================================
% ASSEMBLY THE STIFFNESS MATRIX FOR EACH SPECIES
%===============================================================================

  A1 = -bim1a_advection_diffusion(x, mobility(1)*Vth(1), 1, 1, -valence(1)*phi/Vth(1));
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

  df1_dn1=[];
  df1_dn2=[];
  df1_dn3=[];
  df1_dn4=[];
  df1_dn5=[];
  df1_dn6=[];

  df2_dn1=[];
  df2_dn2=[];
  df2_dn3=[];
  df2_dn4=[];
  df2_dn5=[];
  df2_dn6=[];

  df3_dn1=[];
  df3_dn2=[];
  df3_dn3=[];
  df3_dn4=[];
  df3_dn5=[];
  df3_dn6=[];

  df4_dn1=[];
  df4_dn2=[];
  df4_dn3=[];
  df4_dn4=[];
  df4_dn5=[];
  df4_dn6=[];

  df5_dn1=[];
  df5_dn2=[];
  df5_dn3=[];
  df5_dn4=[];
  df5_dn5=[];
  df5_dn6=[];

  df6_dn1=[];
  df6_dn2=[];
  df6_dn3=[];
  df6_dn4=[];
  df6_dn5=[];
  df6_dn6=[];
  for ii=1:N
    J_N=compute_change_rates_jacobian(n(ii,:),reactions, index);
    df1_dn1=[df1_dn1; J_N(1,1)];
    df1_dn2=[df1_dn2; J_N(1,2)];
    df1_dn3=[df1_dn3; J_N(1,3)];
    df1_dn4=[df1_dn4; J_N(1,4)];
    df1_dn5=[df1_dn5; J_N(1,5)];
    df1_dn6=[df1_dn6; J_N(1,6)];

    df2_dn1=[df2_dn1; J_N(2,1)];
    df2_dn2=[df2_dn2; J_N(2,2)];
    df2_dn3=[df2_dn3; J_N(2,3)];
    df2_dn4=[df2_dn4; J_N(2,4)];
    df2_dn5=[df2_dn5; J_N(2,5)];
    df2_dn6=[df2_dn6; J_N(2,6)];

    df3_dn1=[df3_dn1; J_N(3,1)];
    df3_dn2=[df3_dn2; J_N(3,2)];
    df3_dn3=[df3_dn3; J_N(3,3)];
    df3_dn4=[df3_dn4; J_N(3,4)];
    df3_dn5=[df3_dn5; J_N(3,5)];
    df3_dn6=[df3_dn6; J_N(3,6)];

    df4_dn1=[df4_dn1; J_N(4,1)];
    df4_dn2=[df4_dn2; J_N(4,2)];
    df4_dn3=[df4_dn3; J_N(4,3)];
    df4_dn4=[df4_dn4; J_N(4,4)];
    df4_dn5=[df4_dn5; J_N(4,5)];
    df4_dn6=[df4_dn6; J_N(4,6)];

    df5_dn1=[df5_dn1; J_N(5,1)];
    df5_dn2=[df5_dn2; J_N(5,2)];
    df5_dn3=[df5_dn3; J_N(5,3)];
    df5_dn4=[df5_dn4; J_N(5,4)];
    df5_dn5=[df5_dn5; J_N(5,5)];
    df5_dn6=[df5_dn6; J_N(5,6)];

    df6_dn1=[df6_dn1; J_N(6,1)];
    df6_dn2=[df6_dn2; J_N(6,2)];
    df6_dn3=[df6_dn3; J_N(6,3)];
    df6_dn4=[df6_dn4; J_N(6,4)];
    df6_dn5=[df6_dn5; J_N(6,5)];
    df6_dn6=[df6_dn6; J_N(6,6)];
  endfor

  J=[A1 - diag(Mk*df1_dn1),- diag(Mk*df1_dn2),- diag(Mk*df1_dn3),- diag(Mk*df1_dn4),- diag(Mk*df1_dn5),- diag(Mk*df1_dn6), zeros(N, N);
    - diag(Mk*df2_dn1), A2- diag(Mk*df2_dn2),- diag(Mk*df2_dn3), - diag(Mk*df2_dn4),- diag(Mk*df2_dn5),- diag(Mk*df2_dn6), zeros(N, N);
    - diag(Mk*df3_dn1), - diag(Mk*df3_dn2), A3- diag(Mk*df3_dn3),- diag(Mk*df3_dn4),- diag(Mk*df3_dn5),- diag(Mk*df3_dn6), zeros(N, N);
    - diag(Mk*df4_dn1), - diag(Mk*df4_dn2),- diag(Mk*df4_dn3), A4- diag(Mk*df4_dn4),- diag(Mk*df4_dn5),- diag(Mk*df4_dn6), zeros(N, N);
    - diag(Mk*df5_dn1), - diag(Mk*df5_dn2),- diag(Mk*df5_dn3),- diag(Mk*df5_dn4), A5- diag(Mk*df5_dn5),- diag(Mk*df5_dn6), zeros(N, N);
    - diag(Mk*df6_dn1), - diag(Mk*df6_dn2),- diag(Mk*df6_dn3),- diag(Mk*df6_dn4),- diag(Mk*df6_dn5), A6- diag(Mk*df6_dn6), zeros(N, N);
    -diag(q*Mk*n1* valence(1)),-diag(q*Mk*n2* valence(2)),-diag(q*Mk*n3* valence(3)),-diag(q*Mk*n4* valence(4)),-diag(q*Mk*n5* valence(5)), -diag(q*Mk*n6* valence(6)), P];


  J=sparse(J);
%===============================================================================
% ASSEMBLY THE JACOBIAN (dR_du')
%===============================================================================
  M=zeros((N_species+1)*N,(N_species+1)*N);
  for k=1:N_species
    M(1+(N*(k-1)): k*N, 1+(N*(k-1)):k*N) = Mk;
  endfor
  DJac=sparse(M);

endfunction
