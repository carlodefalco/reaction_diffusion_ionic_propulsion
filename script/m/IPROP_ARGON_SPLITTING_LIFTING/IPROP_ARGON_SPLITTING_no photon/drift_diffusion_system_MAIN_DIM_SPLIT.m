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
[r,idx] = read_reactions(file_in_loadpath ("balcon_et_al_argon_ionization.json"));
%pretty_print_reactions (r);
N_species = numfields(idx);
elements = fieldnames(idx);
%load reactions.gz

%===============================================================================
% INPUT DATA (INTERNATIONAL SYSTEM OF UNITS (SI) ADOPTED)
%===============================================================================

%% MESH GENERATION (1D)
L = 2e-3; %[m]
N = 200;
x = linspace (0, L, N)';

%% TIME EXTREMA FOR INTEGRATION
T0 = 0;%[s]
T1 = 8e-5;%[s];
%T_vec=[0 logspace(-9, log10(T1), 1e3)];
T_vec=linspace(0, T1, 100);

%%con numeri tamnburini simulo 1 millesimo di secondo con 1 milione di punti nel vettore dei tempi
%% con i numeri calcolati per avere pressione di 1 atm simulo ?? con 1 milione di punti nel vettore dei tempi
%% FUNDAMENTAL CONSTANTS
q       = 1.6e-19;      %elementary charge [C]
k_b     = 1.380649e-23; %Boltzmann constant [J K^-1]
epsilon = 8.8e-12;      %dielectric constant in vacuum [C^2 N−1 m−2]
me      = 9.11e-31;      %[Kg]
%%GAS PROPERTIES
Mean_Energy = 4.92; %[eV]
%T_e         = ((2/3)*Mean_Energy*1.60218e-19/k_b);  %[K]
T_e=10000;
T_ions      = 300;    %[K]
pressure    = 101325; %[Pa]

%% MOBILITY AND DIFFUSION COEFFICIENTS
mobility = zeros (numfields (idx)-1, 1);
##mobility(idx.("e"))    = 4.23e-2;%[m^2 / (V*s)]
##mobility(idx.("Ar"))   = 1e-4;   %[m^2 / (V*s)]
##mobility(idx.("Ar+"))  = 1e-4;   %[m^2 / (V*s)]
##mobility(idx.("Ar*"))  = 1e-4;   %[m^2 / (V*s)]
##mobility(idx.("Ar2+")) = 1e-4;   %[m^2 / (V*s

mobility(idx.("e"))    = 1e-1;%[cm^2 / (V*s)]1e-1;%[m^2 / (V*s)]
mobility(idx.("Ar"))   = 1e-3;%[cm^2 / (V*s)]1e-3;%[m^2 / (V*s)]
mobility(idx.("Ar+"))  = 1e-3;%[cm^2 / (V*s)]1e-3;%[m^2 / (V*s)]
mobility(idx.("Ar*"))  = 1e-3;%[cm^2 / (V*s)]1e-3;%[m^2 / (V*s)]
mobility(idx.("Ar2+")) = 1e-3;%[cm^2 / (V*s)]1e-3;%[m^2 / (V*s)]


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
%Vth(idx.("h_nu")) = Vth_ions;

%% APPLIED POTENTIAL
V= 5000; %[V]
%V0 = @(t) V*[zeros(size(t)); ones(size(t))];
V0=V*[0 1];
%% VALENCE NUMBER
valence = zeros (numfields (idx)-1, 1);
valence(idx.("e"))    = -1;              %[1]
valence(idx.("Ar"))   = 0;               %[1]
valence(idx.("Ar+"))  = +1;              %[1]
valence(idx.("Ar*"))  = 0;               %[1]
valence(idx.("Ar2+")) = +1;              %[1]

%===============================================================================
% INITIAL CONDITIONS
%===============================================================================
##%% INITIAL LOCAL DENSITIES
##state0 = zeros (numfields (idx), 1);
##state0(idx.("Ar"))   = 2.5e+25;     %[m^-3]`
##state0(idx.("e"))    = 1.0001e+17;  %[m^-3]
##state0(idx.("Ar+"))  = 1.0e17;      %[m^-3]
##state0(idx.("Ar2+")) = 1.0e13;      %[m^-3]
##state0(idx.("Ar*"))  = 1.0e17;      %[m^-3]
##state0(idx.("h_nu")) = 0;           %[m^-3]
##
##y0 = zeros(N*(N_species),1);
##a=20;
##%distribution=@(x) (1- ((x-.9*L)/(.1*L)).^ 2) .* (x>.8*L & x<L)
##distribution=@(x) exp(-( (x-L).^2 *(a/L^2)))
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = state0(idx.("e")).*distribution(x);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
##y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = state0(idx.("Ar2+")).*distribution(x).* (x>0) ;
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = state0(idx.("Ar+")).*distribution(x).* (x>0);
##y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*")) = state0(idx.("Ar*")).*distribution(x);

%===============================================================================
% TEST 2 (DIR  NON OMOGENEE)
%===============================================================================
%% INITIAL LOCAL DENSITIES
state0 = zeros (numfields (idx), 1);
state0(idx.("Ar"))   = 2.5e+25;%[m^-3]`
state0(idx.("e"))    = 1.0e+18;%[m^-3]
state0(idx.("Ar+"))  = 0;%[m^-3]
state0(idx.("Ar2+")) = 0;%[m^-3]
state0(idx.("Ar*"))  = 0;%[m^-3]
state0(idx.("h_nu")) = 0;%[m^-3]

y0 = zeros(N*(N_species),1);
y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = state0(idx.("e"));
y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
%y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = state0(idx.("Ar+"));

initial_el_approx=zeros(size(x));
for m=1:1000
  initial_el_approx+=-(18/(m*pi))*state0(1)*(-1)^m.*sin(m*pi.*x/L);
endfor
##figure()
##plot(x, initial_el_approx)
##hold on
##plot(x, 9*state0(1)/L.*x)
%===============================================================================
% BOUNDARY CONDITIONS
%===============================================================================
% BC is a matrix of dimensions (N_species+1, 2)to assign 2 boundary values
% at each species and at the potential phi;
##BC_species = zeros(N_species,2);
##BC_species(idx.("Ar"),:)= ones(1,2)*state0(idx.("Ar"));
##BC_species(idx.("e"),:)=y0([1+(idx.("e")-1)*N idx.("e")*N]);
##BC_species(idx.("Ar2+"),:)=y0([1+(idx.("Ar2+")-1)*N idx.("Ar2+")*N]);
##BC_species(idx.("Ar+"),:)=y0([1+(idx.("Ar+")-1)*N idx.("Ar+")*N]);
##BC_species(idx.("Ar*"),:)=y0([1+(idx.("Ar*")-1)*N idx.("Ar*")*N]);
BC_species = zeros(N_species,2);
BC_species(idx.("Ar"),:)= ones(1,2)*state0(idx.("Ar"));
BC_species(idx.("e"),:)=state0(idx.("e")).*[1 10];
BC_species(idx.("Ar2+"),:)=state0(idx.("Ar2+")).*[1 10];
BC_species(idx.("Ar+"),:)=state0(idx.("Ar+")).*[1 10];
BC_species(idx.("Ar*"),:)=state0(idx.("Ar*")).*[1 10];
BC_phi= V0;
nodal_chem=[1 N N+1 2*N 2*N+1 3*N 3*N+1 4*N 4*N+1 5*N 5*N+1 6*N ];
nodal_trans=[1 N N+1 2*N 2*N+1 3*N 3*N+1 4*N 4*N+1 5*N 5*N+1 N_species*N];
##y0=u(1:N*N_species, end);
##T_vec=linspace(t, T1, 10000-columns(u));
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
y0([nodal_chem N_species*N+1 (N_species+1)*N])=[];

%===============================================================================
% INTEGRATION IN TIME
%===============================================================================
##tol = @(t) .001 * (1-((T_vec(1) - t) ./ (T_vec(1) - T_vec(end))).^.4) + 1e-6   * ((T_vec (1) - t) ./ (T_vec(1) - T_vec(end))).^.4;
##tol = @(t) 1e-1;
##mag =V*ones((N_species)*N, 1);
##max_el=max((state0(idx.("e"))).*distribution(x));
##max_Ar_ion=max((state0(idx.("Ar+"))).*distribution(x));
##max_Ar_ex=max((state0(idx.("Ar*"))).*distribution(x));
##max_Ar2_ion=max((state0(idx.("Ar2+"))).*distribution(x));
##
##mag(1:N) = max_el;
##mag(1+N:2*N) = state0(2);
##mag(1+2*N:3*N) = max_Ar_ion;
##mag(1+3*N:4*N) = max_Ar_ex;
##mag(1+4*N:5*N) = max_Ar2_ion;
##mag(1+5*N:6*N) = max_Ar_ex;

tol = @(t) 1e-5;
mag = V*ones((N_species)*(N), 1);

max_el=max(10*state0(idx.("e")));
max_Ar_ion=max(10*state0(idx.("Ar+")));
max_Ar_ex=max_el;
max_Ar2_ion=max_el;

mag(1:N) = max_el;
mag(1+N:2*N) = state0(2);
mag(1+2*N:3*N) = max_el;
mag(1+3*N:4*N) = max_Ar_ex;
mag(1+4*N:5*N) = max_Ar2_ion;
mag(nodal_trans) = [];
M=kron(eye(N_species-1), Mk);  %Mass matrix for all species ((N-2)*N_species x (N-2)*N_species)
mass1=sparse(((N_species)*N),((N_species)*N));
mass1(1:((N_species-1)*N), 1:((N_species-1)*N))=M;
mass2=eye(N_species*N);

M(nodal_trans(1:end-2),:)=[];
M(:, nodal_trans(1:end-2))=[];
mass1(nodal_trans,:)=[];
mass1(:, nodal_trans)=[];
mass2(nodal_chem,:)=[];
mass2(:,nodal_chem)=[];

mass={mass1; mass2};
BC=[BC_species(1,1);BC_species(1,2);BC_species(2,1);BC_species(2,2);BC_species(3,1);BC_species(3,2);BC_species(4,1);BC_species(4,2);BC_species(5,1);BC_species(5,2);BC_phi(1);BC_phi(2) ];
rate={@(y) assembly_rate_transport( y, x, r,idx, valence, mobility,BC,N, q, epsilon, Vth, P, Mk,k_b, me, T_e,nodal_trans);
      @(y) assembly_rate_chem( y, x, r,idx, valence, mobility,BC_species, BC_phi,N, q, epsilon, Vth, P, Mk)};
jacrate={@(y) assembly_jacrate_transport( y, x, r,idx, valence, mobility,BC,N, q, epsilon, Vth, P, Mk, k_b, me, T_e,nodal_trans);
         @(y) assembly_jacrate_chem( y, x, r,idx, valence, mobility,BC_species, BC_phi,N, q, epsilon, Vth, P, Mk)};

warning ('off', 'Octave:nearly-singular-matrix')
[u, error] = integrate_adaptive (y0, T_vec, tol, mag, rate, jacrate, mass);
save -binary -z results.gz

%===============================================================================
% DATA PROCESSING
%===============================================================================
%% rearrenge y vector to plot results by species
y=u';
N=N-2;
n1 =[ones(numel(T_vec),1).*BC(1) y(:,1:N) ones(numel(T_vec),1).*BC(2)];
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
for ii =1:1:numel(t)

 semilogy(x, n1(ii, :),'LineWidth', 1,'r')
## hold on
## semilogy(x, n2(ii, :)*1e-7,'LineWidth', 1,'g')
## semilogy(x, n3(ii, :),'LineWidth', 1,'b')
## semilogy(x, n4(ii, :),'LineWidth', 1,'m')
## semilogy(x, n5(ii, :),'LineWidth', 1,'c')
 %semilogy(x, n6(ii, :),'LineWidth', 1,'y')
 axis ([min(x) max(x) 1e5 2.5e21]);
 str=sprintf('time = %d s', t(ii));
 legend('e','Ar*1e-7', 'Ar+', 'Ar*', 'Ar2+', 'h_nu', 'Location', 'NorthWest')
 title(str)
 hold off
 drawnow
endfor
keyboard
el_exact_staz= (BC_species(1,2)-BC_species(1,1)).*x/L + BC_species(1,1);
ionAr_exact_staz=(BC_species(3,2)-BC_species(3,1)).*x/L + BC_species(3,1);
n_el_exact=zeros(numel(t), numel(x));
n_ar_exact=zeros(numel(t), numel(x));

el_exact_trans= @(t,x,m) +(BC_species(1,2)-BC_species(1,1))*(2/(pi*m)).*((-1)^(m)).*sin(m*pi.*x/L).*exp(-(m*pi/L)^2*mobility(1)*Vth_e*t);

for k=1:numel(t)
  n_el_exact(k,:)=el_exact_staz;
  for m=1:10000
   n_el_exact(k,:)+=el_exact_trans(t(k),x, m)';
  endfor
endfor
err_abs=zeros(size(t));
err_rel=zeros(size(t));
for ii = 1 :1: numel(t)
 err_abs(ii)=norm(n1(ii,:)-n_el_exact(ii,:), inf);
 err_rel(ii)=err_abs(ii)/norm(n_el_exact(ii,:),inf);
endfor
figure()
plot(t, err_abs)
title('err abs')
ylabel('norm(n-nex,inf)')
xlabel('t')
figure()
plot(t, err_rel)
title('err rel, N=200')
ylabel('norm(n-nex,inf)/norm(nex,inf)')
xlabel('t')

phi_ex=zeros(numel(t), numel(x));
C=-q*L/epsilon*(2*state0(1))+5000/L;
phi_ex_staz=q/epsilon*(9*state0(1).*x.^3/(6*L)+state0(1).*x.^2/2)+C*x;
phi_trans= @(t,x,m) -(q/epsilon)*(18*state0(1)*L^2/(m*pi)^3).*((-1)^(m)).*sin(m*pi.*x/L).*exp(-(m*pi/L)^2*mobility(1)*Vth_e*t);
for k=1:numel(t)
  phi_ex(k,:)=phi_ex_staz;
  for m=1:1000
   phi_ex(k,:)+=phi_trans(t(k),x, m)';
  endfor
endfor
figure()
for ii=1:1
  plot(x, phi_ex_staz)
  hold on
  plot(x, y(end,1001:1200), 'x')
endfor

figure()
for ii = 1 :10: 40

 plot(x, n1(ii, :),'LineWidth', 0.5,'rx')
 hold on
 plot(x, n_el_exact(ii, :),'LineWidth', 1,'b')
 %semilogy(x, n3(ii, :),'LineWidth', 1,'b')
 %semilogy(x, n_ar_exact(ii, :),'LineWidth', 1,'bx')

 axis ([min(x) max(x) 1e2 2e19]);
 str=sprintf('time = %d s', t(ii));
 legend('e','ex', 'Ar+', 'Ar_ex')
 title(str)
 hold off
 drawnow
endfor
