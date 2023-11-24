
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
addpath (canonicalize_file_name ("./READ_REACTIONS"));
##[r,idx] = read_reactions(file_in_loadpath ("balcon_et_al_argon_ionization.json"));
[r,idx] = read_reactions(file_in_loadpath ("ionization.json"));
%pretty_print_reactions (r);
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
T1 = 1e-5;%[s];
T_vec=linspace(0, T1, 10000);

%%con numeri tamnburini simulo 1 millesimo di secondo con 1 milione di punti nel vettore dei tempi
%% con i numeri calcolati per avere pressione di 1 atm simulo ?? con 1 milione di punti nel vettore dei tempi
%% FUNDAMENTAL CONSTANTS
q       = 1.6e-19;      %elementary charge [C]
k_b     = 1.380649e-23; %Boltzmann constant [J K^-1]
epsilon = 8.854e-12;    %dielectric constant in vacuum [C^2 N−1 m−2]
me      = 9.11e-31;     %[Kg]
%%GAS PROPERTIES
Mean_Energy = 5.36; %[eV]
T_e         = ((2/3)*Mean_Energy*1.60218e-19/k_b);  %[K]
T_ions      = 300;    %[K]
pressure    = 101325; %[Pa]

%% MOBILITY AND DIFFUSION COEFFICIENTS
##mobility = zeros (numfields (idx)-1, 1);
mobility = zeros (numfields (idx), 1);
mobility(idx.("e"))    = 2.6475e-1;%[m^2 / (V*s)]
mobility(idx.("Ar"))   = 1e-3;   %[m^2 / (V*s)]
mobility(idx.("Ar+"))  = 1e-3;   %[m^2 / (V*s)]
##mobility(idx.("Ar*"))  = 1e-4;   %[m^2 / (V*s)]
##mobility(idx.("Ar2+")) = 1e-4;   %[m^2 / (V*s)]

%% Diffusion coefficient according to Einstein–Smoluchowski
%% relations is D_k= mu_k Vth
%% Vth= k_b*T/q
Vth_e =k_b*T_e/q;
Vth_ions = k_b*T_ions/q;

##Vth=zeros (numfields (idx)-1, 1); %ON if BALCON
Vth=zeros (numfields (idx), 1); %ON If Townsend

Vth(idx.("e"))    = Vth_e;
Vth(idx.("Ar"))   = Vth_ions;%*1e3;
Vth(idx.("Ar+"))  = Vth_ions;%*1e3;
##Vth(idx.("Ar*"))  = Vth_ions;%*1e3;
##Vth(idx.("Ar2+")) = Vth_ions;%*1e3;
%Vth(idx.("h_nu")) = Vth_ions;

%% APPLIED POTENTIAL
V= 500; %[V]
%V0 = @(t) V*[zeros(size(t)); ones(size(t))];
V0=V*[0 1];
%% VALENCE NUMBER
##valence = zeros (numfields (idx)-1, 1);%ON if BALCON
valence = zeros (numfields (idx), 1);%ON If Townsend
valence(idx.("e"))    = -1;              %[1]
valence(idx.("Ar"))   = 0;               %[1]
valence(idx.("Ar+"))  = +1;              %[1]
##valence(idx.("Ar*"))  = 0;               %[1]
##valence(idx.("Ar2+")) = +1;              %[1]

%===============================================================================
% INITIAL CONDITIONS
%===============================================================================
%% INITIAL LOCAL DENSITIES
state0 = zeros (numfields (idx), 1);
state0(idx.("Ar"))   = 2.5e+25;     %[m^-3]`
state0(idx.("e"))    = 4e+18;  %[m^-3]
state0(idx.("Ar+"))  = 4e+18;      %[m^-3]
##state0(idx.("Ar2+")) = 1.0e17;      %[m^-3]
##state0(idx.("Ar*"))  = 1.0e19;      %[m^-3]
##state0(idx.("h_nu")) = 0;           %[m^-3]

y0 = zeros(N*(N_species),1);
A=7;
##E0=2.5e5;
##b=1*(k_b*T_e/(2*pi*me))^0.5+mobility(1)*E0;
##x0=L*(1-b/(2*A*mobility(1)*Vth_e));
##x0=1.99999999e-3;
##A=(b/(2*mobility(1)*Vth_e)) * (L/(L-x0));
%distribution=@(x) (1- ((x-.9*L)/(.1*L)).^ 2) .* (x>.8*L & x<L)
distribution=@(x) exp(-( (x-L).^2 *(A/L^2)))
y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = state0(idx.("e")).*distribution(x) .*(x>0);
y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = state0(idx.("Ar2+")).*distribution(x).* (x>0) ;
y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = state0(idx.("Ar+")).*distribution(x).* (x>0);
##y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*")) = state0(idx.("Ar*")).*distribution(x);

%===============================================================================
% BOUNDARY CONDITIONS
%===============================================================================
% BC is a matrix of dimensions (N_species+1, 2)to assign 2 boundary values
% at each species and at the potential phi;
BC_species = zeros(2,1);
BC_phi= V0;

##nodal_trans=[(idx.("Ar+")-1)*N+1 (idx.("Ar2+")-1)*N+1 5*N+1 N_species*N ];%ON If Balcon
nodal_DIR=[(idx.("Ar+")-1)*N+1 (N_species)*N+1 (N_species+1)*N ];%ON If Townsend
nodal_ROB=[idx.("e")*N];
nodal={nodal_DIR;nodal_ROB};


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
for [value, key]=idx
  if (idx.(key)== "h_nu")
      n_0(:, end)=[]
  endif
endfor
b=n_0*valence;
f_phi_no_bc=q*Mk(2:end-1, 2:end-1)*b;
f_phi = f_phi_no_bc - P(2:end-1,1).*BC_phi(1)  -  P(2:end-1,end).*BC_phi(2);
phi_0(2:end-1) = P(2:end-1, 2:end-1) \ f_phi;
y0=[y0; phi_0];

%===============================================================================
% INTEGRATION IN TIME
%===============================================================================
tol = @(t) 1e-1;
mag =V*ones((N_species+1)*N, 1);
max_el=max((state0(idx.("e"))).*distribution(x));
max_Ar_ion=max((state0(idx.("Ar+"))).*distribution(x));
##max_Ar_ex=max((state0(idx.("Ar*"))).*distribution(x));
##max_Ar2_ion=max((state0(idx.("Ar2+"))).*distribution(x));

mag(1:N) = max_el;
mag(1+N:2*N) = state0(2);
mag(1+2*N:3*N) = max_Ar_ion;
##mag(1+3*N:4*N) = max_Ar_ex;
##mag(1+4*N:5*N) = max_Ar2_ion;
%mag(1+5*N:6*N) = max_Ar_ex;

%%COMPLETE
##M=kron(eye(N_species-1), Mk);  %Mass matrix for all species ((N-2)*N_species x (N-2)*N_species)
##mass1=sparse((N_species)*N,(N_species)*N);
##mass1(1:(N_species-1)*N, 1:(N_species-1)*N)=M;
##mass2=M;

%% townsend simplified
M=kron(eye(N_species), Mk);  %Mass matrix for all species ((N-2)*N_species x (N-2)*N_species)
mass1=sparse((N_species+1)*N,(N_species+1)*N);
mass1(1:(N_species)*N, 1:(N_species)*N)=M;
mass2=M;
mass={mass1; mass2};

##BC=[BC_species(1,1);BC_species(2,1);BC_phi(1);BC_phi(2) ];
BC=[BC_species(1,1);BC_phi(1);BC_phi(2) ];
##boundary_DIR=zeros((numel(y0)-N),1);
##for k=1:numel(nodal{1})
##    boundary_DIR(nodal{1}(k))=BC(k);
##endfor
boundary_ROB=(1*(k_b*T_e/(2*pi*me))^0.5);
keyboard
boundary_values={BC; boundary_ROB};
current=@(n, phi) compute_flux(n, phi,q, Vth_e, mobility(1),x);
rate={@(y) assembly_rate_transport( y, x, r,idx, valence, mobility,q, epsilon, Vth, P, Mk);
      @(y,current) assembly_rate_chem( y, x, r,idx, valence, mobility,N, q, epsilon, Vth, P, Mk, current)};
jacrate={@(y) assembly_jacrate_transport( y, x, r,idx, valence, mobility,q, epsilon, Vth, P, Mk);
         @(y) assembly_jacrate_chem( y, x, r,idx, valence, mobility,N, q, epsilon, Vth, P, Mk)};

warning ('off', 'Octave:nearly-singular-matrix')
[u, error] = integrate_adaptive (y0, T_vec, tol, mag, rate, jacrate, mass, nodal, boundary_values,x, current);
save -binary -z results.gz
keyboard
%===============================================================================
% DATA PROCESSING
%===============================================================================
%% rearrenge y vector to plot results by species

figure
semilogy(x, u(1:200, end-3), '--', 'LineWidth', 1.5)
hold on
semilogy(x, u(1:200, end-2), '--', 'LineWidth', 1.5)
semilogy(x, u(1:200, end-1))
semilogy(x, u(1:200, end))
ylim([1e1 1e18])
legend("first", "second", "third", "last")
title("electron density (last for iteration. final time 5e-4")
print ("-dpng", sprintf ("electron_density.png"))

figure
semilogy(x, u(401:600, end-3), '--', 'LineWidth', 1.5)
hold on
semilogy(x, u(401:600, end-2), '--', 'LineWidth', 1.5)
semilogy(x, u(401:600, end-1))
semilogy(x, u(401:600, end))
legend("first", "second", "third", "last")
title("Ar+ density (last for iteration. final time 5e-4")
print ("-dpng", sprintf ("Ar+_density.png"))

figure
semilogy(x, u(801:1000, end-3),'--', 'LineWidth', 1.5)
hold on
semilogy(x, u(801:1000, end-2),'--', 'LineWidth', 1.5)
semilogy(x, u(801:1000, end-1))
semilogy(x, u(801:1000, end))
legend("first", "second", "third", "last")
title("Ar2+ density (last for iteration. final time 5e-4")
print ("-dpng", sprintf ("Ar2+.png"))

figure
semilogy(x, u(801:1000, end-3),'--', 'LineWidth', 1.5)
hold on
semilogy(x, u(801:1000, end-2),'--', 'LineWidth', 1.5)
semilogy(x, u(801:1000, end-1))
semilogy(x, u(801:1000, end))
legend("first", "second", "third", "last")
title("Ar* density (last for iteration. final time 5e-4")
print ("-dpng", sprintf ("Arex.png"))
figure
plot(x, gradient([0 ;u(1201:end,end-3); 500],x),'--', 'LineWidth', 1.5)
hold on
plot(x, gradient([0 ;u(1201:end,end-2); 500],x),'--', 'LineWidth', 1.5)
plot(x, gradient([0; u(1201:end,end-1); 500],x))
plot(x, gradient([0; u(1201:end,end); 500],x))
legend("first", "second", "third", "last")
title("E(V/m)(last for iteration. final time 5e-4")
print ("-dpng", sprintf ("E_field.png"))
y=u';
n1 = y(:,1:N);
n2 = y(:,N+1:2*N);
n3 = y(:,2*N+1:3*N);
n4 = y(:,3*N+1:4*N) ;
n5 = y(:,4*N+1:5*N) ;
n6 = y(:,5*N+1:6*N) ;
%phi = y(:,6*N+1:7*N) ;
##V = V0(t')(2, :);
##for ii = 1: 100
##  plot (t, u(:, ii))
##  print ("-dpng", sprintf ("frame_%4.4d.png", ii))
##endfor
t=T_vec;
figure()
for ii =1:100:numel(t)

 semilogy(x, n1(ii, :),'LineWidth', 1,'r')
 hold on
 semilogy(x, n2(ii, :)*1e-7,'LineWidth', 1,'g')
 semilogy(x, n3(ii, :),'LineWidth', 1,'b')
 semilogy(x, n4(ii, :),'LineWidth', 1,'m')
 semilogy(x, n5(ii, :),'LineWidth', 1,'c')
 semilogy(x, n6(ii, :),'LineWidth', 1,'y')
 axis ([min(x) max(x) 1e5 2.5e21]);
 str=sprintf('time = %d s', t(ii));
 legend('e','Ar*1e-7', 'Ar+', 'Ar*', 'Ar2+', 'h_nu', 'Location', 'NorthWest')
 title(str)
 hold off
 drawnow
endfor
load results.gz
figure
for ii=1:2:2
semilogy(x, u_chemical(1:200, 7))
hold on
semilogy(x, u_chemical(1:200, 7))
drawnow
endfor
