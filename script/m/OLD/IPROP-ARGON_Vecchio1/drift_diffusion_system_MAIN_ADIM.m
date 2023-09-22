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
pretty_print_reactions (r);;
N_species = numfields(idx);
elements = fieldnames(idx);

%===============================================================================
% INPUT DATA (INTERNATIONAL SYSTEM OF UNITS (SI) ADOPTED)
%===============================================================================

%% MESH GENERATION (1D)
L = 2e-3; %[m]
N = 1000;
x = linspace (0, L, N)';

%% TIME EXTREMA FOR INTEGRATION
T0 = 0;%[s]
T1 = 1.02e-10;%[s]
T_vec=linspace(T0, T1, 1e2);

%% FUNDAMENTAL CONSTANTS
q       = 1.6e-19;      %elementary charge [C]
k_b     = 1.380649e-23; %Boltzmann constant [J K^-1]
epsilon = 8.8e-12;      %dielectric constant in vacuum [C V−1 m−1]

%%GAS PROPERTIES
T_e       = 10000;  %[K]
T_ions    = 300;    %[K]
pressure  = 101325; %[Pa]
m_e=9.31e-31; %[Kg]

%% MOBILITY AND DIFFUSION COEFFICIENTS
mobility = zeros (numfields (idx), 1);
mobility(idx.("e"))    = 1e-2;;%[m^2 / (V*s)]
mobility(idx.("Ar"))   = 1e-6;%[m^2 / (V*s)]
mobility(idx.("Ar+"))  = 1e-6;%[m^2 / (V*s)]
mobility(idx.("Ar*"))  = 1e-6;%[m^2 / (V*s)]
mobility(idx.("Ar2+")) = 1e-6;%[m^2 / (V*s)]

%% Diffusion coefficient according to Einstein–Smoluchowski
%% relations is D_k= mu_k Vth
%% Vth= k_b*T/q
Vth_e = k_b*T_e/q;
Vth_ions = k_b*T_ions/q;
Vth=[Vth_e Vth_ions];

%% APPLIED POTENTIAL
V= 5000; %[V]
V0 = @(t) V*[zeros(size(t)); ((t-T1/100)*1e6 .* ((t>=T1/100)&(t<=T1/100)) + 1.0 .* (t>T1/100))];

%% VALENCE NUMBER
valence = zeros (numfields (idx), 1);
valence(idx.("e"))    = -1;              %[1]
valence(idx.("Ar"))   = 0;               %[1]
valence(idx.("Ar+"))  = +1;              %[1]
valence(idx.("Ar*"))  = 0;               %[1]
valence(idx.("Ar2+")) = +1;              %[1]

%===============================================================================
% INITIAL CONDITIONS
%===============================================================================
%% INITIAL LOCAL DENSITIES
state0 = zeros (numfields (idx), 1);
state0(idx.("Ar"))   = 2.5e+25;
state0(idx.("e"))    = 1.0e+19;
state0(idx.("Ar+"))  = 1.0e+19;
state0(idx.("Ar2+")) = 0;
state0(idx.("Ar*"))  = 0;
state0(idx.("h_nu")) = 0;

y0 = zeros(N*(N_species),1);
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = x .* (L - x) * state0(idx.("e"))/ (L/2)^2;
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = x .* (L - x) * state0(idx.("Ar+"))/ (L/2)^2;
%keyboard
##mu=0.95*L;
##sigma=sqrt(5e-10);
##distribution= @(x) 1/(sqrt(2*pi*sigma^2))*exp(-0.5*((x-mu)./sigma).^2);
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = state0(idx.("e")).*distribution(x);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = state0(idx.("Ar+")).*distribution(x);

distribution=@(x) (1- ((x-.95*L)/(.05*L)).^ 2) .* (x>.9*L)
y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = (state0(idx.("e"))/0.066543).*distribution(x);
y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = (state0(idx.("Ar+"))/0.066543).*distribution(x);

%===============================================================================
% BOUNDARY CONDITIONS
%===============================================================================
% BC is a matrix of dimensions (N_species+1, 2)to assign 2 boundary values
% at each species and at the potential phi;
BC_species = zeros(N_species-1,2);
BC_species(idx.("Ar"),:)= ones(1,2)*state0(idx.("Ar"));
BC_species(idx.("e"),:)=y0([1+(idx.("e")-1)*N idx.("e")*N]);
BC_species(idx.("Ar+"),:)=y0([1+(idx.("Ar+")-1)*N idx.("Ar+")*N]);
BC_phi= V0;


%===============================================================================
% REFERENCE QUANTITIES
%===============================================================================
x_bar = L;

%phi_bar= Vth_e;
phi_bar= V;

v_bar = (k_b*T_e/m_e)^0.5;
%v_bar = max(mobility)*phi_bar/x_bar;

n_bar = state0(idx.("Ar"));



t_bar= x_bar/v_bar;
R_bar=(n_bar*v_bar)/x_bar;

%===============================================================================
% NON DIMENSIONAL QUANTITIES
%===============================================================================
L_hat=L/x_bar;
epsilon_hat= (epsilon*phi_bar)/(x_bar^2*q*n_bar);
state0_hat=state0./n_bar;
x_hat=x./x_bar;
V_hat=V/phi_bar;
Vth_hat=Vth./phi_bar;
mobility_hat=mobility.*(phi_bar/(x_bar*v_bar));
R_bar=((n_bar*v_bar)/x_bar);
T_vec_hat=T_vec./t_bar;


%===============================================================================
% BOUNDARY CONDITIONS(NON DIMENSIONAL)
%===============================================================================
% BC is a matrix of dimensions (N_species+1, 2)to assign 2 boundary values
% at each species and at the potential phi;
BC_species_hat = BC_species./n_bar;
BC_phi_hat= @(t) (V/phi_bar)*[zeros(size(t)); ((t-T1/100)*1e6 .* ((t>=T1/100)&(t<=T1/100)) + 1.0 .* (t>T1/100))];
%===============================================================================
% INITIAL CONDITIONS (NON DIMENSIONAL)
%===============================================================================

y0_hat = y0./n_bar;

%===============================================================================
% COMPUTE CONSISTENT INITIAL CONDITIONS
%===============================================================================
y0dot=consistent_initial_conditions_ADIM(y0_hat,T0,x_hat,Vth_hat,mobility_hat,valence,epsilon_hat, r,idx,BC_species, BC_phi, R_bar, n_bar);



y0_bc=[];
y0dot_bc=[];
for k=1:N_species
  y0_bc    = [y0_bc;y0_hat(2+(k-1)*N: k*N-1)];  %remove left and right boundary from y0
  y0dot_bc = [y0dot_bc;y0dot(2+(k-1)*N: k*N-1)]; %remove left and right boundary from y0dot
endfor
y0dot_hat=y0dot_bc;
y0_hat=y0_bc;

%===============================================================================
% SYSTEM ASSEMBLY
%===============================================================================

fun=@(t,y, ydot) compute_drift_diffusion_reaction_system_ADIM (t, y, ydot, r, idx,...
                       x_hat, N, epsilon_hat, Vth_hat, mobility_hat, valence, BC_species_hat, BC_phi_hat,R_bar,n_bar);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6,...
                'Jacobian', @(t, y, ydot) compute_jacobian_system_ADIM(t, y, ydot,...
              r, idx, x_hat, N, epsilon_hat, Vth_hat, mobility_hat, valence, BC_species_hat, BC_phi_hat, R_bar, n_bar),...
              'MaxOrder', 5);


%===============================================================================
% INTEGRATION IN TIME
%===============================================================================
[t, y] = ode15i (fun, T_vec_hat, y0_hat, y0dot_hat, options);
keyboard
%===============================================================================
% DATA PROCESSING
%===============================================================================
%% rearrenge y vector to plot results by species

n1=zeros(numel(t), N);
n2=zeros(numel(t), N);
n3=zeros(numel(t), N);
n4=zeros(numel(t), N);
n5=zeros(numel(t), N);
n6=zeros(numel(t), N);
phi=zeros(numel(t), N);
for ii=1:numel(t)
  phi_bound=BC_phi(t(ii));
  n1(ii, :)  =[ BC_species(1,1) y(ii,1:(N-2))           BC_species(1,2)];
  n2(ii, :)  =[ BC_species(2,1) y(ii,(N-2)+1:2*(N-2))   BC_species(2,2)];
  n3(ii, :)  =[ BC_species(3,1) y(ii,2*(N-2)+1:3*(N-2)) BC_species(3,2)];
  n4(ii, :)  =[ BC_species(4,1) y(ii,3*(N-2)+1:4*(N-2)) BC_species(4,2)];
  n5(ii, :)  =[ BC_species(5,1) y(ii,4*(N-2)+1:5*(N-2)) BC_species(5,2)];
  n6(ii, :)  =[ BC_species(6,1) y(ii,5*(N-2)+1:6*(N-2)) BC_species(6,2)];
  phi(ii, :) =[ phi_bound(1)    y(ii,6*(N-2)+1:7*(N-2)) phi_bound(2)];

endfor

V = V0(t')(2, :);
%verify conservation of mass
%compute initial total mass
##y0 = [];
##
##for k=1:length(elements)
##  y0=[y0; x.* (L - x)* state0(idx.(elements{k}))/ (L/2)^2];
##endfor
##n10 = y0(1:N);
##n20 = y0(N+1:2*N);
##n30 = y0(2*N+1:3*N);
##n40 = y0(3*N+1:4*N);
##n50 = y0(4*N+1:5*N);
##n60 = y0(5*N+1:6*N);
##
##tot_mass_distribution=sum(n10+n20+n30+n40+n50+n60);
##mass_distribution_computed=sum((n1+n2+n3+n4+n5),2);
##
##for ii=1:100
##  diff=mass_distribution_computed(ii)-tot_mass_distribution
##
##  if diff<1e-7
##    printf("Mass conserved")
##  else
##    printf("Mass not Conserved")
##  endif
##endfor

figure()
for ii = 1 :1e3: numel (t)

  subplot(1, 2, 1)
  %for k= 1:length(elements)
  semilogy(x,n1(ii, 1:N),'LineWidth', 1)
  hold on
  semilogy(x,n2(ii, 1:N),'LineWidth', 1)
  semilogy(x,n3(ii, 1:N),'LineWidth', 1)
  semilogy(x,n4(ii, 1:N),'LineWidth', 1)
  semilogy(x,n5(ii, 1:N),'LineWidth', 1)
  semilogy(x,n6(ii, 1:N),'LineWidth', 1)
  %endfor
  title (sprintf ("%g", t(ii)));
  legend (elements);
  %axis ([min(x) max(x) 0 max(max(state0))]);
  subplot(1, 2, 2)
  plot (t, V, t(ii), V(ii), 'ro')
  title ('V')
  xlabel ('t [s]')
  ylabel ('V [V]')
  drawnow
endfor

##figure()
##semilogy(x,n1(1, 1:N)*2.5e19,'LineWidth', 1,':')
##hold on
##semilogy(x,n2(1, 1:N)*2.5e19,'LineWidth', 1, '-.')
##semilogy(x,n3(1, 1:N)*2.5e19,'LineWidth', 1, '--')
##semilogy(x,n4(1, 1:N)*2.5e19,'LineWidth', 1, '-')
##semilogy(x,n5(1, 1:N)*2.5e19,'LineWidth', 1, '.')
##semilogy(x,n6(1, 1:N)*2.5e19,'LineWidth', 1)
##legend (elements)
##title('initialzation')
##xlabel ('t [s]')
##ylabel ('V [V]')
##drawnow
##
##figure()
##semilogy(x,n1(end, 1:N)*2.5e19,'LineWidth', 1,':')
##hold on
##semilogy(x,n2(end, 1:N)*2.5e19,'LineWidth', 1, '-.')
##semilogy(x,n3(end, 1:N)*2.5e19,'LineWidth', 1, '--')
##semilogy(x,n4(end, 1:N)*2.5e19,'LineWidth', 1,'-')
##semilogy(x,n5(end, 1:N)*2.5e19,'LineWidth', 1,'.')
##semilogy(x,n6(end, 1:N)*2.5e19,'LineWidth', 1)
##legend(elements)
##title('final state')
##xlabel ('t [s]')
##ylabel ('V [V]')
##drawnow


