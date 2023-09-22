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
N = 200;
x = linspace (0, L, N)';

%% TIME EXTREMA FOR INTEGRATION
T0 = 0;%[s]
T1 = 1e-9;%[s];
%T_vec=[0 logspace(-9, log10(T1), 1e3)];
T_vec=linspace(0, T1, 1000);

%%con numeri tamnburini simulo 1 millesimo di secondo con 1 milione di punti nel vettore dei tempi
%% con i numeri calcolati per avere pressione di 1 atm simulo ?? con 1 milione di punti nel vettore dei tempi
%% FUNDAMENTAL CONSTANTS
q       = 1.6e-19;      %elementary charge [C]
k_b     = 1.380649e-23; %Boltzmann constant [J K^-1]
epsilon = 8.8e-12;   %   dielectric constant in vacuum [C^2 N−1 m−2]

%%GAS PROPERTIES
T_e       = 10000;  %[K] (0.8617328149741 eV)
T_ions    = 300;    %[K]
pressure  = 101325; %[Pa]

%% MOBILITY AND DIFFUSION COEFFICIENTS
mobility = zeros (numfields (idx), 1);
mobility(idx.("e"))    = 3.15e-2;%1e-1;%[m^2 / (V*s)]
mobility(idx.("Ar"))   = 1e-4;%[m^2 / (V*s)]
mobility(idx.("Ar+"))  = 1e-4;%[m^2 / (V*s)]
mobility(idx.("Ar*"))  = 1e-4;%[m^2 / (V*s)]
mobility(idx.("Ar2+")) = 1e-4;%[m^2 / (V*s)]

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
%V0 = @(t) V*[zeros(size(t)); ((t-T1/10)*10 .* ((t>=T1/10)&(t<=T1/10)) + 1.0 .* (t>T1/10))];
V0 = @(t) V*[zeros(size(t)); ones(size(t))];

%% VALENCE NUMBER
valence = zeros (numfields (idx), 1);
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
state0(idx.("Ar"))   = 2.5e+25;%[m^-3]`
state0(idx.("e"))    = 1.01e+18;%[m^-3]
state0(idx.("Ar+"))  = 1.0e18;%[m^-3]
state0(idx.("Ar2+")) = 1.0e16;%[m^-3] %PLASMA NON QUASI NEUTRAL--> mettere zero?
state0(idx.("Ar*"))  = 1.0e18;%[m^-3]

y0 = zeros(N*(N_species),1);
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = x .* (L - x) * state0(idx.("e"))/ (L/2)^2,
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = x .* (L - x) * state0(idx.("Ar+"))/ (L/2)^2;
%keyboard
##mu=0.95*L;
##sigma=sqrt(1e-12);
##distribution= @(x) 1/(sqrt(2*pi*sigma^2))*exp(-0.5*((x-mu)./sigma).^2);
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = state0(idx.("e")).*distribution(x);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = state0(idx.("Ar+")).*distribution(x);

distribution=@(x) (1- ((x-.9*L)/(.1*L)).^ 2) .* (x>.8*L & x<L)
y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = (state0(idx.("e"))/0.066543).*distribution(x);
y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = (state0(idx.("Ar2+"))/0.066543).*distribution(x);
y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = (state0(idx.("Ar+"))/0.066543).*distribution(x);
y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*")) = (state0(idx.("Ar*"))/0.066543).*distribution(x);
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = state0(idx.("e")).* (x>.8*L & x<L);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = state0(idx.("Ar2+")).* (x>.8*L & x<L);
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = state0(idx.("Ar+")).* (x>.8*L & x<L);
##y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*")) = state0(idx.("Ar*")).* (x>.8*L & x<L);
%y0(N*idx.("e"))=8.9051e+15;
%y0(N*idx.("Ar+"))=3.4472e+11;
%y0(N*idx.("Ar*"))=2.7835e+15;
%y0(N*idx.("Ar2+"))=8.9048e+15;
##distribution=@(x) (1- ((x-.9*L)/(.1*L)).^ 2) .* (x>.8*L & x<L)
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = (state0(idx.("e"))/0.066543).*distribution(x) + 2e11*exp(1e4*(x-0.8*L)).*(x<=0.8*L& x>0);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##%y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = (state0(idx.("Ar+"))/0.066543).*distribution(x)+ x .* (L - x) * 1e6/ (L/2)^2.*(x<=.8*L & x>0);
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = (state0(idx.("Ar+"))/0.066543).*distribution(x) + 2e11*exp(1e4*(x-0.8*L)).*(x<=0.8*L& x>0);

%===============================================================================
% BOUNDARY CONDITIONS
%===============================================================================
% BC is a matrix of dimensions (N_species+1, 2)to assign 2 boundary values
% at each species and at the potential phi;
BC_species = zeros(N_species,2);
BC_species(idx.("Ar"),:)= ones(1,2)*state0(idx.("Ar"));
BC_species(idx.("e"),:)=y0([1+(idx.("e")-1)*N idx.("e")*N]);
BC_species(idx.("Ar2+"),:)=y0([1+(idx.("Ar2+")-1)*N idx.("Ar2+")*N]);
BC_species(idx.("Ar+"),:)=y0([1+(idx.("Ar+")-1)*N idx.("Ar+")*N]);
BC_species(idx.("Ar*"),:)=y0([1+(idx.("Ar*")-1)*N idx.("Ar*")*N]);

##BC_species(idx.("Ar"),:)= ones(1,2)*state0(idx.("Ar"));
##BC_species(idx.("e"),:)=[y0(1+(idx.("e")-1)*N) 8.9051e+15];
##BC_species(idx.("Ar2+"),:)=[y0(1+(idx.("Ar2+")-1)*N) 8.9048e+15];
##BC_species(idx.("Ar+"),:)=[y0(1+(idx.("Ar+")-1)*N) 3.4472e+11];
##BC_species(idx.("Ar*"),:)=[y0(1+(idx.("Ar*")-1)*N) 2.7835e+15];


BC_phi= V0;
keyboard
%===============================================================================
% COMPUTE CONSISTENT INITIAL CONDITIONS
%===============================================================================
##y0dot=consistent_initial_conditions(y0,T0,x,Vth,mobility,valence,epsilon,q,r,idx,BC_species, BC_phi);
##
##y0_bc=[];
##y0dot_bc=[];
##for k=1:N_species
##  y0_bc    = [y0_bc;y0(2+(k-1)*N: k*N-1)];  %remove left and right boundary from y0
##  y0dot_bc = [y0dot_bc;y0dot(2+(k-1)*N: k*N-1)]; %remove left and right boundary from y0dot
##endfor
##res=compute_drift_diffusion_reaction_system_3(T0, y0_bc, y0dot_bc, r, idx,...
##                       x, N, q, epsilon, Vth, mobility, valence, BC_species, BC_phi);
%===============================================================================
% COMPUTE PHI(0)
%===============================================================================
Mk=bim1a_reaction (x, 1, 1); %Mass matrix for each species (NxN)
P= bim1a_laplacian (x, epsilon, 1); %stiffness matrix for phi(NxN)
phi_0 = zeros (N, 1);
BC_phi0=BC_phi(0);
phi_0([1 end])= BC_phi0;
y0_bc=y0;
y0_bc([1:N:(N_species)*N, N:N:(N_species)*N])=[];
n_0=reshape(y0_bc, N-2, N_species);
b=n_0*valence;
f_phi_no_bc=q*Mk(2:end-1, 2:end-1)*b;
f_phi = f_phi_no_bc - P(2:end-1,1).*BC_phi0(1)  -  P(2:end-1,end).*BC_phi0(2);
phi_0(2:end-1) = P(2:end-1, 2:end-1) \ f_phi;
y0=[y0; phi_0];
y0_bc=[y0_bc; phi_0(2:end-1)];
keyboard

%===============================================================================
% INTEGRATION IN TIME
%===============================================================================



%[t, y] = ode15i (fun,T_vec, y0_bc, y0dot_bc, options);
%[t,y, res_abs, res_rel]=time_integrator_for_drift(T_vec, x, y0, r, idx, valence, mobility,BC_species, BC_phi, N, q, epsilon, Vth, P, Mk);
%tol = @(t) .001 * (1-((T_vec(1) - t) / (T_vec(1) - T_vec(end)))^.4) + 1e-6   * ((T_vec (1) - t) / (T_vec(1) - T_vec(end)))^.4;
tol = @(t) 1e-4;
mag = 5000*ones((N_species+1)*(N-2), 1);
max_el=max((state0(idx.("e"))/0.066543).*distribution(x));
max_Ar_ion=max((state0(idx.("Ar+"))/0.066543).*distribution(x));
max_Ar_ex=max((state0(idx.("Ar*"))/0.066543).*distribution(x));
max_Ar2_ion=max((state0(idx.("Ar2+"))/0.066543).*distribution(x));

mag(1:(N-2)) = max_el;
mag(1+(N-2):2*(N-2)) = state0(2);
mag(1+2*(N-2):3*(N-2)) = max_Ar_ion;
mag(1+3*(N-2):4*(N-2)) = max_Ar_ex;
mag(1+4*(N-2):5*(N-2)) = max_Ar2_ion;
mag(1+5*(N-2):6*(N-2)) = max_el;

M=kron(eye(N_species), Mk(2:end-1, 2:end-1));  %Mass matrix for all species ((N-2)*N_species x (N-2)*N_species)
mass=sparse((N_species+1)*(N-2),(N_species+1)*(N-2));
mass(1:N_species*(N-2), 1:N_species*(N-2))=M;
rate=@(y,t) assembly_rate(t, y, x, r,idx, valence, mobility,BC_species, BC_phi,N, q, epsilon, Vth, P, Mk);
jacrate=@(y,t) assembly_jacrate(t, y, x, r,idx, valence, mobility,BC_species, BC_phi,N, q, epsilon, Vth, P, Mk);
keyboard
warning ('off', 'Octave:nearly-singular-matrix')
u = integrate_adaptive (y0_bc, T_vec, tol, mag, rate, jacrate, mass);
keyboard
%===============================================================================
% DATA PROCESSING
%===============================================================================
%% rearrenge y vector to plot results by species

y=u';
n1 =[ones(size(T_vec'))*BC_species(1,1) y(:,1:N-2)             ones(size(T_vec'))*BC_species(1,2)];
n2 =[ones(size(T_vec'))*BC_species(2,1) y(:,(N-2)+1:2*(N-2))   ones(size(T_vec'))*BC_species(2,2)];
n3 =[ones(size(T_vec'))*BC_species(3,1) y(:,2*(N-2)+1:3*(N-2)) ones(size(T_vec'))*BC_species(3,2)];
n4 =[ones(size(T_vec'))*BC_species(4,1) y(:,3*(N-2)+1:4*(N-2)) ones(size(T_vec'))*BC_species(4,2)];
n5 =[ones(size(T_vec'))*BC_species(5,1) y(:,4*(N-2)+1:5*(N-2)) ones(size(T_vec'))*BC_species(5,2)];
n6 =[zeros(size(T_vec'))*BC_species(6,1) y(:,5*(N-2)+1:6*(N-2)) zeros(size(T_vec'))*BC_species(6,2)];
phi =[zeros(size(T_vec')) y(:,6*(N-2)+1:7*(N-2)) ones(size(T_vec'))*5000];
V = V0(t')(2, :);

t=T_vec;
figure()
for ii = 1 :50: numel(t)

 semilogy(x, n1(ii, :),'LineWidth', 1,'r')
 hold on
 semilogy(x, n2(ii, :)*1e-5,'LineWidth', 1,'b')
 semilogy(x, n3(ii, :),'LineWidth', 1,'b')
 semilogy(x, n4(ii, :),'LineWidth', 1,'m')
 semilogy(x, n5(ii, :),'LineWidth', 1,'c')
  semilogy(x, n6(ii, :),'LineWidth', 1,'y')
 axis ([min(x) max(x) 1e2 2.5e20]);
 str=sprintf('time = %d s', t(ii));
 legend('e','Ar*1e-12', 'Ar+', 'Ar*', 'Ar2+', 'h_nu')
 title(str)
 hold off
 drawnow
endfor

figure()
for ii = 1 :50: numel(t)

 plot(x, phi(ii, :),'LineWidth', 1,'r')
 axis ([min(x) max(x) 0 5030]);
 str=sprintf('time = %d s', t(ii));
 legend('phi')
 title(str)
 hold off
 drawnow
endfor


##figure()
##loglog(t(2:end), (res_rel),'LineWidth', 1,'c')
##axis ([t(2) t(end) 1e-21 2.5e-5]);
##str=sprintf('relative residue');
##legend( ' res Ar' )
##title(str)

##electron_flux=@(x)(+mobility(1)*Vth_e*state0(1)*2*(x-0.9*L)/(0.066543*0.1^2*L^2)...
##                    +mobility(1)*(1-((x-0.9*L)/(0.1*L)).^2)*state0(1)*V/(L*0.066543)).*(x>0.8*L);
##div_electron=@(x)(+mobility(1)*Vth_e*state0(1)*2/(0.066543*0.1^2*L^2)...
##                    -mobility(1)*2*(x-0.9*L)/(0.066543*0.1^2*L^2)*state0(1)*V/L).*(x>0.8*L);
##Fe=electron_flux(x);
##divergence=div_electron(x);
##s_e0=y0(N*(idx.("e")-1)  +1:N*idx.("e")) ;
##s_Ar0=y0(N*(idx.("Ar")-1) +1:N*idx.("Ar")) ;
##s_Ar_excited0=y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*"));
##s_Ar_plus0=y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) ;
##s_Ar2_plus0=y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+"));
##s_hnu0=y0(N*(idx.("h_nu")-1)+1:N*idx.("h_nu"));
##R1_f0 = r(1).rate_coeffs(1).*s_e0.*s_Ar0;
##R2_f0 = r(2).rate_coeffs(1).*s_e0.*s_Ar0;
##R3_f0 = r(3).rate_coeffs(1).*s_e0.*s_Ar_excited0;
##R4_f0 = r(4).rate_coeffs(1).*s_Ar_excited0.^2;
##R5_f0 = r(5).rate_coeffs(1).*s_Ar_plus0.*s_Ar0.^2;
##R6_f0 = r(6).rate_coeffs(1).*s_e0.*s_Ar2_plus0;
##R7_f0 = r(7).rate_coeffs(1).*s_Ar_excited0;
##R8_f0 = r(8).rate_coeffs(1).*s_e0.*s_Ar0;
##s_e0_dot = R1_f0 + R3_f0 + R4_f0 - R6_f0;
##


