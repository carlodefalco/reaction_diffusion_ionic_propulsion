%% script test per verificare che la rimozione del ciclo for in compute_change_rates sia corretta

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
state0 = zeros (numfields (idx), 1);
state0(idx.("Ar"))   = 2.5e+25;%[m^-3]`
state0(idx.("e"))    = 1.01e+18;%[m^-3]
state0(idx.("Ar+"))  = 1.0e18;%[m^-3]
state0(idx.("Ar2+")) = 1.0e16;%[m^-3] %PLASMA NON QUASI NEUTRAL--> mettere zero?
state0(idx.("Ar*"))  = 1.0e12;%[m^-3]
y0 = zeros(N*(N_species),1);
distribution=@(x) (1- ((x-.9*L)/(.1*L)).^ 2) .* (x>.8*L & x<L)
y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = (state0(idx.("e"))/0.066543).*distribution(x);
y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = (state0(idx.("Ar2+"))/0.066543).*distribution(x);
y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = (state0(idx.("Ar+"))/0.066543).*distribution(x);
y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*")) = (state0(idx.("Ar*"))/0.066543).*distribution(x);

%tic
dstate1=compute_change_rates2(y0, r, idx, x);
%toc
%keyboard
%tic
dstate2=compute_change_rates3(y0, r, idx, x);
%toc

keyboard
tic
jac1=compute_change_rates_jacobian2(y0, r, idx, x);
toc
keyboard
tic
jac2=compute_change_rates_jacobian3(y0, r, idx, x);
toc
printf(num2str(sum(sum(dstate1-dstate2))))
printf(num2str(sum(sum(jac1-jac2))))
