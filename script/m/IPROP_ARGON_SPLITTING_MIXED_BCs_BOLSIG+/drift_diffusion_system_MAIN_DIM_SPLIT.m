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
close all
addpath (canonicalize_file_name ("../../../data"));
addpath (canonicalize_file_name ("./Assembly"));
addpath (canonicalize_file_name ("./Integrate"));
addpath (canonicalize_file_name ("./Chemical_Equations"));
addpath (canonicalize_file_name ("./READ_REACTIONS"));
addpath (canonicalize_file_name ("../BolsigPlus"));
%===============================================================================
% LOAD FPL MSH BIM PACKAGEs FOR FE MATRIX ASSEMBLY
%===============================================================================
pkg load fpl bim msh
%===============================================================================
% LOAD DATABASE FOR COEFFICIENTS COMPUTATION
%===============================================================================
load("database.txt")

%% MESH GENERATION (1D)
L = 2e-3; %[m]
N = 200;
x = linspace (0, L, N)';

%% TIME EXTREMA FOR INTEGRATION
T0 = 0;%[s]
T1 = 5e-4;%[s];
%T_vec=[0 logspace(-9, log10(T1), 1e3)];
T_vec=linspace(0, T1, 10000);

%===============================================================================
% READ THE REACTION FILE INPUT
%===============================================================================
load reactions.gz
N_species = numfields(idx);
elements = fieldnames(idx);
for ii=1:numel(r)
 guess=r(ii).rate_coeffs;
 r(ii).rate_coeffs=[];
 r(ii).rate_coeffs=[ones(N,1).*guess(1) ones(N,1).*guess(2)];
endfor
%===============================================================================
% INPUT DATA (INTERNATIONAL SYSTEM OF UNITS (SI) ADOPTED)
%===============================================================================
%%con numeri tamnburini simulo 1 millesimo di secondo con 1 milione di punti nel vettore dei tempi
%% con i numeri calcolati per avere pressione di 1 atm simulo ?? con 1 milione di punti nel vettore dei tempi
%% FUNDAMENTAL CONSTANTS
q       = 1.6e-19;      %elementary charge [C]
k_b     = 1.380649e-23; %Boltzmann constant [J K^-1]
epsilon = 8.8e-12;      %dielectric constant in vacuum [C^2 N−1 m−2]
me      = 9.11e-31;     %[Kg] mass of the electron

%%GAS PROPERTIES
Mean_Energy=5.2;
T_e         = ((2/3)*Mean_Energy*1.60218e-19/k_b);  %[K]
T_ions    = 300;    %[K]
pressure  = 101325; %[Pa]

%% MOBILITY AND DIFFUSION COEFFICIENTS
mobility = zeros (N, numfields (idx)-1);
mobility(:,idx.("e"))    = 4.23e-2*ones(N,1);%[m^2 / (V*s)]
mobility(:,idx.("Ar"))   = 1e-4*ones(N,1);   %[m^2 / (V*s)]
mobility(:,idx.("Ar+"))  = 1e-4*ones(N,1);   %[m^2 / (V*s)]
mobility(:,idx.("Ar*"))  = 1e-4*ones(N,1);   %[m^2 / (V*s)]
mobility(:,idx.("Ar2+")) = 1e-4*ones(N,1);   %[m^2 / (V*s)]
mobility = 2 ./ (1./mobility(1:end-1,:) + 1./mobility(2:end,:));

%% Diffusion coefficient according to Einstein–Smoluchowski
%% relations is D_k= mu_k Vth
%% Vth= k_b*T/q
Vth_e =k_b*T_e/q;
Vth_ions = k_b*T_ions/q;
Vth=zeros (numfields (idx)-1, 1);
Vth(idx.("e"))    = Vth_e;
Vth(idx.("Ar"))   = Vth_ions;
Vth(idx.("Ar+"))  = Vth_ions;
Vth(idx.("Ar*"))  = Vth_ions;
Vth(idx.("Ar2+")) = Vth_ions;
Vth(idx.("h_nu")) = Vth_ions;

%% APPLIED POTENTIAL
V= 500; %[V]
%V0 = @(t) V*[zeros(size(t)); ((t-T1/10)*10 .* ((t>=T1/10)&(t<=T1/10)) + 1.0 .* (t>T1/10))];
%V0 = @(t) V*[zeros(size(t)); ones(size(t))];
V0=V*[0 1];
%% VALENCE NUMBER
valence = zeros (numfields (idx)-1, 1);
valence(idx.("e"))    = -1;              %[1]
valence(idx.("Ar"))   = 0;               %[1]
valence(idx.("Ar+"))  = +1;              %[1]
valence(idx.("Ar*"))  = 0;               %[1]
valence(idx.("Ar2+")) = +1;               %[1]


%===============================================================================
% INITIAL CONDITIONS
%===============================================================================
%% INITIAL LOCAL DENSITIES
state0 = zeros (numfields (idx), 1);
state0(idx.("Ar"))   = 2.5e+25;%[cm^-3]`
state0(idx.("e"))    = 1.05e+17;%[cm^-3]
state0(idx.("Ar+"))  = 1*state0(idx.("e")) ;%[cm^-3]
state0(idx.("Ar2+")) = 0*state0(idx.("e")) ;%[cm^-3]
state0(idx.("Ar*"))  = 0;%[cm^-3]
state0(idx.("h_nu")) = 0;%[cm^-3]

##y0 = zeros(N*(N_species),1);
##a=4;
##distribution=@(x) exp(-( (x-L).^2 *(a/L^2)))
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = state0(idx.("e")).*distribution(x);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = state0(idx.("Ar2+")).*distribution(x).* (x>0) ;
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = state0(idx.("Ar+")).*distribution(x).* (x>0);
##y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*")) = state0(idx.("Ar*")).*distribution(x);
load results_aumento_crush.gz
y0=u(1:1200, end);


##distribution=@(x) (1- ((x-.9*L)/(.1*L)).^ 2) .* (x>.8*L & x<L)
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = (state0(idx.("e"))/0.066543).*distribution(x);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = (state0(idx.("Ar2+"))/0.066543).*distribution(x);
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = (state0(idx.("Ar+"))/0.066543).*distribution(x);
##y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*")) = (state0(idx.("Ar*"))/0.066543).*distribution(x);

%===============================================================================
% BOUNDARY CONDITIONS
%===============================================================================
% BC is a matrix of dimensions (N_species+1, 2)to assign 2 boundary values
% at each species and at the potential phi;
BC_species = zeros(N_species,2);
##BC_species(idx.("Ar"),:)= ones(1,2)*state0(idx.("Ar"));
##BC_species(idx.("e"),:)=y0([1+(idx.("e")-1)*N idx.("e")*N]);
##BC_species(idx.("Ar2+"),:)=y0([1+(idx.("Ar2+")-1)*N idx.("Ar2+")*N]);
##BC_species(idx.("Ar+"),:)=y0([1+(idx.("Ar+")-1)*N idx.("Ar+")*N]);
##BC_species(idx.("Ar*"),:)=y0([1+(idx.("Ar*")-1)*N idx.("Ar*")*N]);
BC_phi= V0;
nodal_trans=[(idx.("Ar+")-1)*N+1 (idx.("Ar2+")-1)*N+1 5*N+1 N_species*N ];
%===============================================================================
% COMPUTE PHI(0)
%===============================================================================
Mk=bim1a_reaction (x, 1, 1); %Mass matrix for each species (NxN)
P= bim1a_laplacian (x, epsilon, 1); %stiffness matrix for phi(NxN)
phi_0 = zeros (N, 1);
phi_0([1 end])= BC_phi;
y0_bc=y0;
y0_bc([1:N:(N_species)*N, N:N:(N_species)*N])=[];
n_0=reshape(y0_bc, N-2, N_species);
n_0(:, end)=[];
b=n_0*valence;
f_phi_no_bc=q*Mk(2:end-1, 2:end-1)*b;
f_phi = f_phi_no_bc - P(2:end-1,1).*BC_phi(1)  -  P(2:end-1,end).*BC_phi(2);
phi_0(2:end-1) = P(2:end-1, 2:end-1) \ f_phi;
y0=[y0; phi_0];


%===============================================================================
% INTEGRATION IN TIME
%===============================================================================

%tol = @(t) .001 * (1-((T_vec(1) - t) / (T_vec(1) - T_vec(end)))^.4) + 1e-6   * ((T_vec (1) - t) / (T_vec(1) - T_vec(end)))^.4;
tol = @(t) 1e-1;
mag = V*ones((N_species)*N, 1);
max_el=max(y0(1:N));
max_Ar=max(y0(N+1:2*N));
max_Ar_ion=max(y0(2*N+1:3*N));
max_Ar_ex=max(y0(3*N+1:4*N));
max_Ar2_ion=max(y0(4*N+1:5*N));

mag(1:N) = max_el;
mag(1+N:2*N) = max_Ar;
mag(1+2*N:3*N) = max_Ar_ion;
mag(1+3*N:4*N) = max_Ar_ex;
mag(1+4*N:5*N) = max_Ar2_ion;



M=kron(eye(N_species-1), Mk);  %Mass matrix for all species (N*N_species x N*N_species)
mass1=sparse((N_species)*N,(N_species)*N);
mass1(1:(N_species-1)*N, 1:(N_species-1)*N)=M;
mass2=eye(N*N_species);

mass={mass1; mass2};
BC=[BC_species(idx.("Ar+"),1);BC_species(idx.("Ar2+"),1);BC_phi(1);BC_phi(2) ];
boundary_values=zeros((numel(y0)-N),1);
for k=1:numel(nodal_trans)
boundary_values(nodal_trans(k))=BC(k);
endfor
Te=@(y) compute_temperature( y,x,idx,N_species, database);
electron_mobility=@(y) compute_electron_mobility( y,x,idx, N_species,database);
rate={@(y) assembly_rate_transport( y, x,idx, valence, mobility,electron_mobility,BC_species, BC_phi,N,N_species, q, epsilon, Vth, P, Mk,me, k_b,Te);
      @(y,r_new) assembly_rate_chem( y, x, r_new,idx, valence, mobility,BC_species, BC_phi,N,N_species, q, epsilon, Vth, P, Mk)};
jacrate={@(y) assembly_jacrate_transport( y, x,idx, valence, mobility,electron_mobility,BC_species, BC_phi,N,N_species, q, epsilon, Vth, P, Mk,me, k_b,Te);
         @(y,r_new) assembly_jacrate_chem( y, x, r_new,idx, valence, mobility,BC_species, BC_phi,N,N_species, q, epsilon, Vth, P, Mk)};

warning ('off', 'Octave:nearly-singular-matrix')
u = integrate_adaptive (y0, T_vec, tol, mag, rate, jacrate, mass,r,idx,x,N_species,database,nodal_trans,boundary_values);
save -binary -z results.gz
%===============================================================================
% DATA PROCESSING
%===============================================================================
%% rearrenge y vector to plot results by species
##phibc=V0(T_vec);
##y=u';
##n1 = y(:,1:N);
##n2 = y(:,N+1:2*N);
##n3 = y(:,2*N+1:3*N);
##n4 = y(:,3*N+1:4*N) ;
##n5 = y(:,4*N+1:5*N) ;
##n6 = y(:,5*N+1:6*N) ;
##phi = y(:,6*N+1:7*N) ;
##V = V0(t')(2, :);
##
##t=T_vec;
##
##figure()
##for ii =1:10:numel(t)
##
## semilogy(x, n1(ii, :),'LineWidth', 1,'r')
## hold on
#### semilogy(x, n2(ii, :)*1e-7,'LineWidth', 1,'b')
## semilogy(x, n3(ii, :),'LineWidth', 1,'b')
#### semilogy(x, n4(ii, :),'LineWidth', 1,'m')
## semilogy(x, n5(ii, :),'LineWidth', 1,'c')
#### semilogy(x, n6(ii, :),'LineWidth', 1,'y')
## axis ([min(x) max(x) 1e10 2.5e21]);
## str=sprintf('time = %d s', t(ii));
## legend('e','Ar*1e-12', 'Ar+', 'Ar*', 'Ar2+', 'h_nu', 'Location', 'NorthWest')
## title(str)
## hold off
## drawnow
##endfor


