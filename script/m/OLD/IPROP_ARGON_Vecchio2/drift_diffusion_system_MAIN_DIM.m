## Copyright (C) 2023 Carlo de Falco
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% Solve :
%%   n_k' = div (mun Vth (grad n_k - grad phi/Vth n_k))+ S_k
%%   0  = -div (epsilon grad (phi)) - q sum_k { z_k  n_k}
%%   n_k(0) = n_k(L)= 0
%%   phi(0) = 0 phi(L) = V0(t)

clc
clear

addpath (canonicalize_file_name ("../../../data"));
addpath (canonicalize_file_name ("./Assembly"));
addpath (canonicalize_file_name ("./Chemical_Equations"));
addpath (canonicalize_file_name ("./READ_REACTIONS"));

%===============================================================================
% LOAD BIM PACKAGE FOR FE MATRIX ASSEMBLY
%===============================================================================
pkg load fpl bim msh

%===============================================================================
% READ THE REACTION FILE INPUT
%===============================================================================
[r,idx] = read_reactions(file_in_loadpath ("balcon_et_al_argon_ionization.json"));
pretty_print_reactions (r);
N_species = numfields(idx);
elements = fieldnames(idx);
%===============================================================================
% INPUT DATA (INTERNATIONAL SYSTEM OF UNITS (SI) ADOPTED)
%===============================================================================

%% MESH GENERATION (1D)
L = 2e-3; %[m]
N =81;
x = linspace (0, L, N)';

%% TIME EXTREMA FOR INTEGRATION
T0 = 0;%[s]
T1 = 100.0e-9;%[s]
T_vec=linspace(T0, T1, 1e2);
T_vec=[0 logspace(-15, -9, 1e2)];
dt=(T1-T0)/numel(T_vec);
%T_vec=[0 logspace(-15, log10 (T1), 1e3)];

%%con numeri tamnburini simulo 1 millesimo di secondo con 1 milione di punti nel vettore dei tempi
%% con i numeri calcolati per avere pressione di 1 atm simulo ?? con 1 milione di punti nel vettore dei tempi
%% FUNDAMENTAL CONSTANTS
q       = 1.6e-19;      %elementary charge [C]
k_b     = 1.380649e-23; %Boltzmann constant [J K^-1]
epsilon = 8.8e-12;      %dielectric constant in vacuum [C V−1 m−1]

%%GAS PROPERTIES
T_e       = 10000;  %[K]
T_ions    = 300;    %[K]
pressure  = 101325; %[Pa]

%% MOBILITY AND DIFFUSION COEFFICIENTS
mobility = zeros (numfields (idx), 1);
mobility(idx.("e"))    = 1e-2;%[m^2 / (V*s)]
mobility(idx.("Ar"))   = 1e-6;%[m^2 / (V*s)]
mobility(idx.("Ar+"))  = 1e-6;%[m^2 / (V*s)]
mobility(idx.("Ar*"))  = 1e-6;%[m^2 / (V*s)]
mobility(idx.("Ar2+")) = 1e-6;%[m^2 / (V*s)]

%% Diffusion coefficient according to Einstein–Smoluchowski
%% relations is D_k= mu_k Vth
%% Vth= k_b*T/q
Vth_e =k_b*T_e/q;
Vth_ions = k_b*T_ions/q;
Vth=zeros (numfields (idx), 1);
Vth(idx.("e"))    = Vth_e;
Vth(idx.("Ar"))   = Vth_ions;
Vth(idx.("Ar+"))  = Vth_ions;
Vth(idx.("Ar*"))  = Vth_ions;
Vth(idx.("Ar2+")) = Vth_ions;
Vth(idx.("h_nu")) = Vth_ions;

%% APPLIED POTENTIAL
V= 5000; %[V]
V0 = @(t) V*[zeros(size(t)); ((t-T1/100)*1e6 .* ((t>=T1/100)&(t<=T1/100)) + 1.0 .* (t>T1/100))];

%% VALENCE NUMBER
valence = zeros (numfields (idx), 1);
valence(idx.("e"))    = -1;              %[1]
valence(idx.("Ar"))   = 0;               %[1]
valence(idx.("Ar+"))  = +1;              %[1]
valence(idx.("Ar*"))  = 0;               %[1]
valence(idx.("Ar2+")) = 1;               %[1]

%===============================================================================
% INITIAL CONDITIONS
%===============================================================================
%% INITIAL LOCAL DENSITIES
state0 = zeros (numfields (idx), 1);
state0(idx.("Ar"))   = 2.5e+25;
state0(idx.("e"))    = 1.0e+12;
state0(idx.("Ar+"))  = 1.0e12;
state0(idx.("Ar2+")) = 0;
state0(idx.("Ar*"))  = 0;
state0(idx.("h_nu")) = 0;

y0 = zeros(N*(N_species),1);
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = x .* (L - x) * state0(idx.("e"))/ (L/2)^2,
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = x .* (L - x) * state0(idx.("Ar+"))/ (L/2)^2;
%keyboard
##mu=0.95*L;
##sigma=sqrt(5e-12);
##distribution= @(x) 1/(sqrt(2*pi*sigma^2))*exp(-0.5*((x-mu)./sigma).^2);
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = state0(idx.("e")).*distribution(x);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = state0(idx.("Ar+")).*distribution(x);

distribution=@(x) (1- ((x-.95*L)/(.05*L)).^ 2) .* (x>.9*L)
y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = (state0(idx.("e"))/0.066543).*distribution(x);
y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = (state0(idx.("Ar+"))/0.066543).*distribution(x);

##L_cm=L*1e2;
##E=V/L_cm;
##p=pressure*0.00750062;
##control=E/p;
##A=12;
##B=180;
##alpha=p*A*exp(-B/(E/p))*1e2;
##ne0=1e6;
##ne= @(x) y0(1:N).*exp(alpha*x);
%===============================================================================
% BOUNDARY CONDITIONS
%===============================================================================
% BC is a matrix of dimensions (N_species+1, 2)to assign 2 boundary values
% at each species and at the potential phi;
BC_species = zeros(N_species,2);
BC_species(idx.("Ar"),:)= ones(1,2)*state0(idx.("Ar"));
BC_species(idx.("e"),:)=y0([1+(idx.("e")-1)*N idx.("e")*N]);
BC_species(idx.("Ar+"),:)=y0([1+(idx.("Ar+")-1)*N idx.("Ar+")*N]);
BC_phi= V0;

##  n0=[];
##  for k=1:N_species
##    n0=[ n0, y0(2+N*(k-1):N*k-1)];
##  endfor
##  chemistry=[];
##  for ii=1:N-2
##    chemistry=[chemistry; compute_change_rates(n0(ii,:),r,idx)] ;
##  endfor
##  chemistry2 =compute_change_rates2(n0, r, idx,x);
##  keyboard
##J=compute_change_rates_jacobian2(n0, r, idx,x);
##  df1_dn1=[];
##  df1_dn2=[];
##  df1_dn3=[];
##  df1_dn4=[];
##  df1_dn5=[];
##  df1_dn6=[];
##
##  df2_dn1=[];
##  df2_dn2=[];
##  df2_dn3=[];
##  df2_dn4=[];
##  df2_dn5=[];
##  df2_dn6=[];
##
##  df3_dn1=[];
##  df3_dn2=[];
##  df3_dn3=[];
##  df3_dn4=[];
##  df3_dn5=[];
##  df3_dn6=[];
##
##  df4_dn1=[];
##  df4_dn2=[];
##  df4_dn3=[];
##  df4_dn4=[];
##  df4_dn5=[];
##  df4_dn6=[];
##
##  df5_dn1=[];
##  df5_dn2=[];
##  df5_dn3=[];
##  df5_dn4=[];
##  df5_dn5=[];
##  df5_dn6=[];
##
##  df6_dn1=[];
##  df6_dn2=[];
##  df6_dn3=[];
##  df6_dn4=[];
##  df6_dn5=[];
##  df6_dn6=[];
##%===============================================================================
##% ASSEMBLY THE JACOBIAN (dR_du)
##%===============================================================================
##  for ii=1:N
##    J_N=compute_change_rates_jacobian( n0(ii,:),r, idx);
##
##
##
##    df1_dn1=[df1_dn1; J_N(1,1)];
##    df1_dn2=[df1_dn2; J_N(1,2)];
##    df1_dn3=[df1_dn3; J_N(1,3)];
##    df1_dn4=[df1_dn4; J_N(1,4)];
##    df1_dn5=[df1_dn5; J_N(1,5)];
##    df1_dn6=[df1_dn6; J_N(1,6)];
##
##    df2_dn1=[df2_dn1; J_N(2,1)];
##    df2_dn2=[df2_dn2; J_N(2,2)];
##    df2_dn3=[df2_dn3; J_N(2,3)];
##    df2_dn4=[df2_dn4; J_N(2,4)];
##    df2_dn5=[df2_dn5; J_N(2,5)];
##    df2_dn6=[df2_dn6; J_N(2,6)];
##
##    df3_dn1=[df3_dn1; J_N(3,1)];
##    df3_dn2=[df3_dn2; J_N(3,2)];
##    df3_dn3=[df3_dn3; J_N(3,3)];
##    df3_dn4=[df3_dn4; J_N(3,4)];
##    df3_dn5=[df3_dn5; J_N(3,5)];
##    df3_dn6=[df3_dn6; J_N(3,6)];
##
##    df4_dn1=[df4_dn1; J_N(4,1)];
##    df4_dn2=[df4_dn2; J_N(4,2)];
##    df4_dn3=[df4_dn3; J_N(4,3)];
##    df4_dn4=[df4_dn4; J_N(4,4)];
##    df4_dn5=[df4_dn5; J_N(4,5)];
##    df4_dn6=[df4_dn6; J_N(4,6)];
##
##    df5_dn1=[df5_dn1; J_N(5,1)];
##    df5_dn2=[df5_dn2; J_N(5,2)];
##    df5_dn3=[df5_dn3; J_N(5,3)];
##    df5_dn4=[df5_dn4; J_N(5,4)];
##    df5_dn5=[df5_dn5; J_N(5,5)];
##    df5_dn6=[df5_dn6; J_N(5,6)];
##
##    df6_dn1=[df6_dn1; J_N(6,1)];
##    df6_dn2=[df6_dn2; J_N(6,2)];
##    df6_dn3=[df6_dn3; J_N(6,3)];
##    df6_dn4=[df6_dn4; J_N(6,4)];
##    df6_dn5=[df6_dn5; J_N(6,5)];
##    df6_dn6=[df6_dn6; J_N(6,6)];
##  endfor
##
##  J2=-[-diag(df1_dn1),-diag(df1_dn2),-diag(df1_dn3),-diag(df1_dn4),-diag(df1_dn5),-diag(df1_dn6);
##      -diag(df2_dn1), -diag(df2_dn2),-diag(df2_dn3), -diag(df2_dn4),-diag(df2_dn5),-diag(df2_dn6);
##      -diag(df3_dn1), -diag(df3_dn2),-diag(df3_dn3),-diag(df3_dn4),-diag(df3_dn5),-diag(df3_dn6);
##      -diag(df4_dn1), -diag(df4_dn2),-diag(df4_dn3),-diag(df4_dn4),-diag(df4_dn5),-diag(df4_dn6);
##      -diag(df5_dn1), -diag(df5_dn2),-diag(df5_dn3),-diag(df5_dn4), -diag(df5_dn5),-diag(df5_dn6);
##      -diag(df6_dn1), -diag(df6_dn2),-diag(df6_dn3),-diag(df6_dn4),-diag(df6_dn5), -diag(df6_dn6);
##    ];

%===============================================================================
% COMPUTE CONSISTENT INITIAL CONDITIONS
%===============================================================================
[y0,y0dot]=consistent_initial_conditions(y0,T0,x,Vth,mobility,valence,epsilon,q,r,idx,BC_species, BC_phi);
keyboard
res=compute_drift_diffusion_reaction_system(T0, y0, y0dot, r, idx,...
                       x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi);
aa=[0
9.9893e+19
1.9727e+20
2.9213e+20
3.8445e+20
4.7424e+20
5.615e+20
6.4624e+20
7.2844e+20
8.0812e+20
8.8526e+20
9.5988e+20
1.032e+21
1.1015e+21
1.1686e+21
1.2331e+21
1.295e+21
1.3545e+21
1.4114e+21
1.4658e+21
1.5176e+21
1.5669e+21
1.6137e+21
1.658e+21
1.6997e+21
1.7389e+21
1.7756e+21
1.8098e+21
1.8414e+21
1.8705e+21
1.897e+21
1.921e+21
1.9425e+21
1.9615e+21
1.978e+21
1.9919e+21
2.0032e+21
2.0121e+21
2.0184e+21
2.0222e+21
2.0235e+21
2.0222e+21
2.0184e+21
2.0121e+21
2.0032e+21
1.9919e+21
1.978e+21
1.9615e+21
1.9425e+21
1.921e+21
1.897e+21
1.8705e+21
1.8414e+21
1.8098e+21
1.7756e+21
1.7389e+21
1.6997e+21
1.658e+21
1.6137e+21
1.5669e+21
1.5176e+21
1.4658e+21
1.4114e+21
1.3545e+21
1.295e+21
1.2331e+21
1.1686e+21
1.1015e+21
1.032e+21
9.5988e+20
8.8526e+20
8.0812e+20
7.2844e+20
6.4624e+20
5.615e+20
4.7424e+20
3.8445e+20
2.9213e+20
1.9727e+20
9.9893e+19
0
-2.071e+27
-7.7148e+30
-1.5234e+31
-2.2559e+31
-2.9688e+31
-3.6621e+31
-4.3359e+31
-4.9902e+31
-5.625e+31
-6.2402e+31
-6.8359e+31
-7.4121e+31
-7.9688e+31
-8.5059e+31
-9.0234e+31
-9.5215e+31
-1e+32
-1.0459e+32
-1.0898e+32
-1.1318e+32
-1.1719e+32
-1.21e+32
-1.2461e+32
-1.2803e+32
-1.3125e+32
-1.3428e+32
-1.3711e+32
-1.3975e+32
-1.4219e+32
-1.4443e+32
-1.4648e+32
-1.4834e+32
-1.5e+32
-1.5146e+32
-1.5273e+32
-1.5381e+32
-1.5469e+32
-1.5537e+32
-1.5586e+32
-1.5615e+32
-1.5625e+32
-1.5615e+32
-1.5586e+32
-1.5537e+32
-1.5469e+32
-1.5381e+32
-1.5273e+32
-1.5146e+32
-1.5e+32
-1.4834e+32
-1.4648e+32
-1.4443e+32
-1.4219e+32
-1.3975e+32
-1.3711e+32
-1.3428e+32
-1.3125e+32
-1.2803e+32
-1.2461e+32
-1.21e+32
-1.1719e+32
-1.1318e+32
-1.0898e+32
-1.0459e+32
-1e+32
-9.5215e+31
-9.0234e+31
-8.5059e+31
-7.9688e+31
-7.4121e+31
-6.8359e+31
-6.2402e+31
-5.625e+31
-4.9902e+31
-4.3359e+31
-3.6621e+31
-2.9688e+31
-2.2559e+31
-1.5234e+31
-7.7148e+30
-2.071e+27
0
-7.7148e+30
-1.5234e+31
-2.2559e+31
-2.9687e+31
-3.6621e+31
-4.3359e+31
-4.9902e+31
-5.625e+31
-6.2402e+31
-6.8359e+31
-7.4121e+31
-7.9687e+31
-8.5059e+31
-9.0234e+31
-9.5215e+31
-1e+32
-1.0459e+32
-1.0898e+32
-1.1318e+32
-1.1719e+32
-1.21e+32
-1.2461e+32
-1.2803e+32
-1.3125e+32
-1.3428e+32
-1.3711e+32
-1.3975e+32
-1.4219e+32
-1.4443e+32
-1.4648e+32
-1.4834e+32
-1.5e+32
-1.5146e+32
-1.5273e+32
-1.5381e+32
-1.5469e+32
-1.5537e+32
-1.5586e+32
-1.5615e+32
-1.5625e+32
-1.5615e+32
-1.5586e+32
-1.5537e+32
-1.5469e+32
-1.5381e+32
-1.5273e+32
-1.5146e+32
-1.5e+32
-1.4834e+32
-1.4648e+32
-1.4443e+32
-1.4219e+32
-1.3975e+32
-1.3711e+32
-1.3428e+32
-1.3125e+32
-1.2803e+32
-1.2461e+32
-1.21e+32
-1.1719e+32
-1.1318e+32
-1.0898e+32
-1.0459e+32
-1e+32
-9.5215e+31
-9.0234e+31
-8.5059e+31
-7.9687e+31
-7.4121e+31
-6.8359e+31
-6.2402e+31
-5.625e+31
-4.9902e+31
-4.3359e+31
-3.6621e+31
-2.9687e+31
-2.2559e+31
-1.5234e+31
-7.7148e+30
0
0
6.4508e+20
1.2738e+21
1.8863e+21
2.4823e+21
3.0621e+21
3.6255e+21
4.1726e+21
4.7034e+21
5.2178e+21
5.7159e+21
6.1977e+21
6.6632e+21
7.1123e+21
7.545e+21
7.9615e+21
8.3616e+21
8.7454e+21
9.1128e+21
9.464e+21
9.7987e+21
1.0117e+22
1.0419e+22
1.0705e+22
1.0975e+22
1.1228e+22
1.1465e+22
1.1685e+22
1.1889e+22
1.2077e+22
1.2248e+22
1.2404e+22
1.2542e+22
1.2665e+22
1.2771e+22
1.2861e+22
1.2934e+22
1.2992e+22
1.3032e+22
1.3057e+22
1.3065e+22
1.3057e+22
1.3032e+22
1.2992e+22
1.2934e+22
1.2861e+22
1.2771e+22
1.2665e+22
1.2542e+22
1.2404e+22
1.2248e+22
1.2077e+22
1.1889e+22
1.1685e+22
1.1465e+22
1.1228e+22
1.0975e+22
1.0705e+22
1.0419e+22
1.0117e+22
9.7987e+21
9.464e+21
9.1128e+21
8.7454e+21
8.3616e+21
7.9615e+21
7.545e+21
7.1123e+21
6.6631e+21
6.1977e+21
5.7159e+21
5.2178e+21
4.7034e+21
4.1726e+21
3.6255e+21
3.0621e+21
2.4824e+21
1.8863e+21
1.2738e+21
6.4508e+20
0
0
7.7148e+30
1.5234e+31
2.2559e+31
2.9688e+31
3.6621e+31
4.3359e+31
4.9902e+31
5.625e+31
6.2402e+31
6.8359e+31
7.4121e+31
7.9688e+31
8.5059e+31
9.0234e+31
9.5215e+31
1e+32
1.0459e+32
1.0898e+32
1.1318e+32
1.1719e+32
1.21e+32
1.2461e+32
1.2803e+32
1.3125e+32
1.3428e+32
1.3711e+32
1.3975e+32
1.4219e+32
1.4443e+32
1.4648e+32
1.4834e+32
1.5e+32
1.5146e+32
1.5273e+32
1.5381e+32
1.5469e+32
1.5537e+32
1.5586e+32
1.5615e+32
1.5625e+32
1.5615e+32
1.5586e+32
1.5537e+32
1.5469e+32
1.5381e+32
1.5273e+32
1.5146e+32
1.5e+32
1.4834e+32
1.4648e+32
1.4443e+32
1.4219e+32
1.3975e+32
1.3711e+32
1.3428e+32
1.3125e+32
1.2803e+32
1.2461e+32
1.21e+32
1.1719e+32
1.1318e+32
1.0898e+32
1.0459e+32
1e+32
9.5215e+31
9.0234e+31
8.5059e+31
7.9687e+31
7.4121e+31
6.8359e+31
6.2402e+31
5.625e+31
4.9902e+31
4.3359e+31
3.6621e+31
2.9688e+31
2.2559e+31
1.5234e+31
7.7148e+30
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0];
                       keyboard
%===============================================================================
% SYSTEM ASSEMBLY
%===============================================================================
fun=@(t,y, ydot) compute_drift_diffusion_reaction_system(t, y, ydot, r, idx,...
                       x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi);
%===============================================================================
% JACOBIAN TEST
%===============================================================================
[J0, DJ0]= compute_jacobian_system(0, y0_bc, y0dot_bc,...
          r, idx, x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi);
keyboard
##function [J2,DJ] = jacobian2 (t, y, ydot, r, idx, x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi)
##
##  T=t;
##  f2=@(y) compute_drift_diffusion_reaction_system_3(T, y, ydot, r, idx,...
##                       x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi);
##  J2=complex_step_diff(y, f2);
##
##  M=bim1a_reaction (x, 1, 1);
##  DJ=sparse(6*(N-2),6*(N-2));
##  DJ(1:N-2, 1:N-2) = M(2:end-1, 2:end-1);
##  DJ((N-2)+1:2*(N-2), (N-2)+1:2*(N-2)) = M(2:end-1, 2:end-1);
##  DJ(2*(N-2)+1:3*(N-2), 2*(N-2)+1:3*(N-2)) = M(2:end-1, 2:end-1);
##  DJ(3*(N-2)+1:4*(N-2), 3*(N-2)+1:4*(N-2)) = M(2:end-1, 2:end-1);
##  DJ(4*(N-2)+1:5*(N-2), 4*(N-2)+1:5*(N-2)) = M(2:end-1, 2:end-1);
##  DJ(5*(N-2)+1:6*(N-2), 5*(N-2)+1:6*(N-2)) = M(2:end-1, 2:end-1);
##
##endfunction
##[J2,DJ2] = jacobian2 (0, y0_bc, y0dot_bc, r, idx, x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi);
##keyboard
%===============================================================================
% INTEGRATION IN TIME
%===============================================================================
##S=diag(ones(N-2,1))+diag(ones(N-3,1),1)+diag(ones(N-3,1),-1);
##A=diag(ones(N-2,1));
##DJPattern=[S A A A A A;
##           A S A A A A;
##           A A S A A A;
##           A A A S A A;
##           A A A A S A;
##           A A A A A A];
##DJPattern=sparse(DJPattern);
##DDJpattern=sparse(eye(N_species*(N-2)));
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'InitialStep',1e-30,...
          'Jacobian',  @(t, y, ydot) compute_jacobian_system(t, y, ydot,...
           r, idx, x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi),...
         'MaxOrder', 5, 'Stats', 'on');
##options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'InitialStep',1e-30,...
##          'JPattern', [DJPattern, DDJpattern],...
##         'MaxOrder', 5, 'Stats', 'on');

[t, y] = ode15i (fun,T_vec, y0_bc, y0dot_bc, options);

keyboard
%===============================================================================
% DATA PROCESSING
%===============================================================================
%% rearrenge y vector to plot results by species
n1 =[ones(size(t))*BC_species(1,1) y(:,1:N-2)             ones(size(t))*BC_species(1,2)];
n2 =[ones(size(t))*BC_species(2,1) y(:,(N-2)+1:2*(N-2))   ones(size(t))*BC_species(2,2)];
n3 =[ones(size(t))*BC_species(3,1) y(:,2*(N-2)+1:3*(N-2)) ones(size(t))*BC_species(3,2)];
n4 =[ones(size(t))*BC_species(4,1) y(:,3*(N-2)+1:4*(N-2)) ones(size(t))*BC_species(4,2)];
n5 =[ones(size(t))*BC_species(5,1) y(:,4*(N-2)+1:5*(N-2)) ones(size(t))*BC_species(5,2)];
n6 =[zeros(size(t)) y(:,5*(N-2)+1:6*(N-2)) zeros(size(t))];

V = V0(t')(2, :);
dt=numel(t)/(1e1);
figure()
for ii = 1 :dt: numel (t)

  subplot(1, 2, 1)
  %for k= 1:length(elements)
  semilogy(x,n1(ii, 1:N),'LineWidth', 1,'r')
  hold on
  semilogy(x,n2(ii, 1:N),'LineWidth', 1,'g')
  semilogy(x,n3(ii, 1:N),'LineWidth', 1,'y')
  semilogy(x,n4(ii, 1:N),'LineWidth', 1,'b')
  semilogy(x,n5(ii, 1:N),'LineWidth', 1,'m')
  semilogy(x,n6(ii, 1:N),'LineWidth', 1,'c')
  %endfor
  hold off
  title (sprintf ("%g", t(ii)));
  legend (elements);
  axis ([min(x) max(x) 10e-10 max(max(state0))]);
  subplot(1, 2, 2)
  plot (t, V, t(ii), V(ii), 'ro')
  title ('V')
  xlabel ('t [s]')
  ylabel ('V [V]')
  drawnow
endfor


figure()
plot(x,n1(1, 1:N),'LineWidth', 1,':')
hold on
plot(x,n3(1, 1:N),'LineWidth', 1, '--')


legend ('e', 'Ar+')
title('initialzation')
xlabel ('t [s]')
ylabel ('V [V]')
drawnow

figure()
semilogy(x,n1(end, 1:N),'LineWidth', 1,':')
hold on
semilogy(x,n2(end, 1:N),'LineWidth', 1, '-.')
semilogy(x,n3(end, 1:N),'LineWidth', 1, '--')
semilogy(x,n4(end, 1:N),'LineWidth', 1,'-')
semilogy(x,n5(end, 1:N),'LineWidth', 1,'.')
semilogy(x,n6(end, 1:N),'LineWidth', 1)
legend(elements)
title('final state')
xlabel ('t [s]')
ylabel ('V [V]')
drawnow


