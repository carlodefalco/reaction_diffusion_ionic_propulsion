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
N = 1000;
x = linspace (0, L, N)';

%% TIME EXTREMA FOR INTEGRATION
T0 = 0;%[s]
T1 = 1.0e-8;%[s];
T_vec=[0 logspace(-10, log10(T1), 5e2)];
%T_vec=[0 logspace(-12, log10 (T1), 1e3)];
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
state0(idx.("Ar+"))  = 1.0e+12;
state0(idx.("Ar2+")) = 0;
state0(idx.("Ar*"))  = 0;
state0(idx.("h_nu")) = 0;

y0 = zeros(N*(N_species),1);
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = x .* (L - x) * state0(idx.("e"))/ (L/2)^2,
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = x .* (L - x) * state0(idx.("Ar+"))/ (L/2)^2;
%keyboard
mu=0.95*L;
sigma=sqrt(1e-12);
distribution= @(x) 1/(sqrt(2*pi*sigma^2))*exp(-0.5*((x-mu)./sigma).^2);
y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = state0(idx.("e")).*distribution(x);
y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = state0(idx.("Ar+")).*distribution(x);

##distribution=@(x) (1- ((x-.95*L)/(.05*L)).^ 2) .* (x>.9*L)
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = (state0(idx.("e"))/0.066543).*distribution(x);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = (state0(idx.("Ar+"))/0.066543).*distribution(x);


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

%===============================================================================
% COMPUTE CONSISTENT INITIAL CONDITIONS
%===============================================================================
y0dot=consistent_initial_conditions(y0,T0,x,Vth,mobility,valence,epsilon,q,r,idx,BC_species, BC_phi);

y0_bc=[];
y0dot_bc=[];
for k=1:N_species
  y0_bc    = [y0_bc;y0(2+(k-1)*N: k*N-1)];  %remove left and right boundary from y0
  y0dot_bc = [y0dot_bc;y0dot(2+(k-1)*N: k*N-1)]; %remove left and right boundary from y0dot
endfor
res=compute_drift_diffusion_reaction_system_3(T0, y0_bc, y0dot_bc, r, idx,...
                       x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi);

%===============================================================================
% SYSTEM ASSEMBLY
%===============================================================================
fun=@(t,y, ydot) compute_drift_diffusion_reaction_system_3(t, y, ydot, r, idx,...
                       x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi);

%===============================================================================
% INTEGRATION IN TIME
%===============================================================================

options = odeset('Jacobian',  @(t, y, ydot) compute_jacobian_system(t, y, ydot,...
           r, idx, x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi),...
         'MaxOrder', 1, 'Stats', 'on');

[t, y] = ode15i (fun,T_vec, y0_bc, y0dot_bc, options);

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

figure()
for ii = 1 :1: numel (t)

 plot(x,n1(ii, :)','LineWidth', 1,'r')
  hold on
  %plot(x,n2(ii, :)','LineWidth', 1,'g')
  plot(x,n3(ii, :)','LineWidth', 1,'y')
  plot(x,n4(ii, :)','LineWidth', 1,'b')
  plot(x,n5(ii, :)','LineWidth', 1,'m')
  plot(x,n6(ii, :)','LineWidth', 1,'c')
  legend ('e', 'Ar+', 'Ar*','Ar2+','h_nu');
  %axis ([min(x) max(x) 1e-5 2.5e25]);
  xlim([0 L])

   hold off
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


