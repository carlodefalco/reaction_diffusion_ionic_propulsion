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
pkg load symbolic
%===============================================================================
% READ THE REACTION FILE INPUT
%===============================================================================
idx.p=1;
idx.n=2;

N_species = numfields(idx);
elements = fieldnames(idx);
%===============================================================================
% INPUT DATA (INTERNATIONAL SYSTEM OF UNITS (SI) ADOPTED)
%===============================================================================

%% MESH GENERATION (1D)
L = 1.5e-6; %[m]
N = 300;
x = linspace (0, L, N)';

%% TIME EXTREMA FOR INTEGRATION
T0 = 0;%[s]
T1 = 2e-10;%[s];
%T_vec=[0 logspace(-9, log10(T1), 1e3)];
T_vec=linspace(0, T1, 1e1);

%%con numeri tamnburini simulo 1 millesimo di secondo con 1 milione di punti nel vettore dei tempi
%% con i numeri calcolati per avere pressione di 1 atm simulo ?? con 1 milione di punti nel vettore dei tempi
%% FUNDAMENTAL CONSTANTS
q       = 1.6e-19;      %elementary charge [C]
epsilon = 4*8.85e-12;   %   dielectric constant in vacuum [C^2 N−1 m−2]

%% MOBILITY AND DIFFUSION COEFFICIENTS
mobility = zeros (numfields (idx), 1);
mobility(idx.("p")) = 1e-1; %[m^2 / (V*s)]
mobility(idx.("n")) = 1e-1; %[m^2 / (V*s)]


Vth=zeros (numfields (idx), 1);
Vth(idx.("p"))   = 26e-3;
Vth(idx.("n"))   = 26e-3;


%% APPLIED POTENTIAL
V= 5000; %[V]
V0 = @(t) V*[zeros(size(t)); zeros(size(t))];
%V0 = @(t) V*[zeros(size(t)); ones(size(t))];

%% VALENCE NUMBER
valence = zeros (numfields (idx), 1);
valence(idx.("p"))   = 1;              %[1]
valence(idx.("n"))   = -1;               %[1]


%===============================================================================
% INITIAL CONDITIONS
%===============================================================================
coeff=1e22;
ui=1e16; %[m^-3]
D= coeff.*(x<.5*L);
A= coeff.*(x>=.5*L);
% INITIAL LOCAL DENSITIES
y0 = zeros(N*(N_species),1);
y0(N*(idx.("p")-1)  +1:N*idx.("p")) = (A(end)+(-A(end)/2+sqrt((A(end)/2)^2+ui^2))) .* (x >  .55*L)+(-D(1)/2+sqrt((D(1)/2)^2+ui^2))*(exp(-.1/Vth(idx.("p")))) .* (x <= .55*L);
y0(N*(idx.("n")-1) +1:N*idx.("n"))  = (D(1)+(-D(1)/2+sqrt((D(1)/2)^2+ui^2))*(exp(-.1/Vth(idx.("n"))))) .* (x <= .45*L)+(-A(end)/2+sqrt((A(end)/2)^2+ui^2)) .* (x >  .45*L);

y0(end)=ui;
y0(1)=ui;

tau=1e-1;
%===============================================================================
% BOUNDARY CONDITIONS
%===============================================================================
% BC is a matrix of dimensions (N_species+1, 2)to assign 2 boundary values
% at each species and at the potential phi;

BC_species = zeros(N_species,2);
%p0=(-D(1)/2 +sqrt(D(1)^2+4*ui^2)/2);
%n0=ui^2/p0;

BC_species(idx.("p"),:)= [ui, y0(N)];
BC_species(idx.("n"),:)= [y0(N+1), ui];

BC_phi= V0;



%===============================================================================
% COMPUTE PHI(0)
%===============================================================================
Mk=bim1a_reaction (x, 1, 1); %Mass matrix for each species (NxN)
P= bim1a_laplacian (x, epsilon, 1); %stiffness matrix for phi(NxN)


%===============================================================================
% REMOVE BOUNDARY VALUES FROM THE INITIAL STATE
%===============================================================================
y0_bc=y0;
y0_bc([1:N:(N_species)*N, N:N:(N_species)*N])=[];
phi_0=zeros(size(x));
y0_bc=[y0_bc; phi_0(2:end-1)];

%===============================================================================
% INTEGRATION IN TIME
%===============================================================================
##[t,y, res_abs, res_rel]=time_integrator_for_drift(T_vec, x, y0, r, idx, valence, mobility,BC_species, BC_phi, N, q, epsilon, Vth, P, Mk);
##keyboard

tol = @(t) .001 * (1-((T_vec(1) - t) / (T_vec(1) - T_vec(end)))^.4) + 1e-6   * ((T_vec (1) - t) / (T_vec(1) - T_vec(end)))^.4;
mag = 100*ones(3*(N-2), 1);
mag(1:(N-2)) = 1e22;
mag((N-2)+1:(N-2)*2) = 1e22;

M=kron(eye(N_species), Mk(2:end-1, 2:end-1));  %Mass matrix for all species ((N-2)*N_species x (N-2)*N_species)
mass=sparse((N_species+1)*(N-2),(N_species+1)*(N-2));
mass(1:N_species*(N-2), 1:N_species*(N-2))=M;
rate=@(y,t) assembly_rate(t, y, x, idx, valence, mobility,BC_species, BC_phi,N, q, epsilon, Vth, P, Mk,ui, tau);
jacrate=@(y,t) assembly_jacrate(t, y, x, idx, valence, mobility,BC_species, BC_phi,N, q, epsilon, Vth, P, Mk,ui, tau);

warning ('off', 'Octave:nearly-singular-matrix')
u_bc = integrate_adaptive (y0_bc, T_vec, tol, mag, rate, jacrate, mass);

u=[ones(size(T_vec)).*BC_species(1,1);
  u_bc(1:N-2,:);
  ones(size(T_vec)).*BC_species(1,2);
  ones(size(T_vec)).*BC_species(2,1);
  u_bc(1+N-2:2*(N-2),:);
  ones(size(T_vec)).*BC_species(2,2);
  zeros(size(T_vec));
  u_bc(2*(N-2)+1:3*(N-2),:);
  zeros(size(T_vec))
  ];
figure()
semilogy(x, u(1:N, end))
hold on
semilogy(x, u(N+1:2*N, end))
legend('p','n')
title('MY CODE')
figure()
plot(x, u(2*N+1:end, :))
title('MY CODE')
