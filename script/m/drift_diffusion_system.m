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
%%   n_k' = div (mun Vth (grad n_k - grad phi/Vth n_k))
%%   0  = -div (epsilon grad (phi)) - q sum_k { z_k  n_k}
%%   n_k(0) = n_k(L)= 0
%%   phi(0) = 0 phi(L) = V0(t)

pkg load fpl bim msh

addpath (canonicalize_file_name ("../../data"));
[r, idx] = read_reactions (file_in_loadpath ("balcon_et_al_argon_ionization.json"));

%pretty_print_reactions (r);

%% initialzation of density, valence number and mobility
state0 = zeros (numfields (idx), 1);
state0(idx.("e"))    = 1.0e+6/2.5e19;
state0(idx.("Ar"))   = 2.5e+19/2.5e19;
state0(idx.("Ar+"))  = 1.0e+6/2.5e19;
state0(idx.("Ar*"))  = 1.0e+10/2.5e19;
state0(idx.("Ar2+")) = 1.0e+3/2.5e19;


%%NOTE THAT THE MOBILITY SHOULD BE COMPUTED NOT IMPOSED: this is for a try
mobility = zeros (numfields (idx), 1);
mobility(idx.("e"))    = 1e-2;
mobility(idx.("Ar"))   = 1e-3;
mobility(idx.("Ar+"))  = 1e-3;
mobility(idx.("Ar*"))  = 1e-3;
mobility(idx.("Ar2+")) = 1e-3;

valence = zeros (numfields (idx), 1);
valence(idx.("e"))    = -1;
valence(idx.("Ar"))   = 0;
valence(idx.("Ar+"))  = 1;
valence(idx.("Ar*"))  = 0;
valence(idx.("Ar2+")) = +1;


elements= fieldnames(idx);
N_species=numfields(idx);
%% Problem Data
L = 2e-3;
N = 81;
x = linspace (0, L, N)';
T0 = 0.99e-4;
T1 = 1.02e-4;
T_vec=linspace(T0, T1, 100);

V0 = @(t) 5000*[zeros(size(t)); ((t-1e-4)*1e6 .* ((t>=1e-4)&(t<=1.01e-4)) + 1.0 .* (t>1.01e-4))];
q  = 1.6e-19; %elementary charge
epsilon = 8.8e-12; %dielectric constant in vacuum
% Diffusion coefficient according to Einstein–Smoluchowski relations is D_k= mu_k Vth
Vth = 26e-3;
%%grid initialzation
y0=[];
for k=1:length(elements)
  y0=[y0; x.* (L - x)* state0(idx.(elements{k}))/ (L/2)^2];
endfor
y0=[y0; zeros(N,1)];
y0dot=zeros(N*N_species+N,1);

%% build the system

fun=@(t,y, ydot) compute_drift_diffusion_reaction_system (t, y, ydot, r, idx, x, N, V0, q, epsilon, Vth, mobility, valence);
options = odeset('RelTol', 10.0^(-4), 'AbsTol', 10.0^(-4), 'InitialStep', 1e-4-T0,...%'normcontrol', 'on', 'refine', 5,...
                'Jacobian', @(t, y, ydot) compute_jacobian_system(t, y, ydot, r, idx, x, N, V0, q, epsilon, Vth, mobility, valence));

%integration in time
[t, y] = ode15i (fun, linspace(T0,T1,100), y0, y0dot, options);

%% rearrenge y vector to plot results by species
n1 = y(:,1:N);
n2 = y(:,N+1:2*N);
n3 = y(:,2*N+1:3*N);
n4 = y(:,3*N+1:4*N);
n5 = y(:,4*N+1:5*N);
n6 = y(:,5*N+1:6*N);
phi= y(:,6*N+1:7*N);
V = V0(t')(2, :);

##figure()
##for ii = 1 : numel (t)
##
##  subplot(1, 2, 1)
##  %for k= 1:length(elements)
##  semilogy(x,n1(ii, 1:N),'LineWidth', 0.2)
##  hold on
##  semilogy(x,n2(ii, 1:N),'LineWidth', 0.2)
##  semilogy(x,n3(ii, 1:N),'LineWidth', 0.2)
##  semilogy(x,n4(ii, 1:N),'LineWidth', 0.2)
##  semilogy(x,n5(ii, 1:N),'LineWidth', 0.2)
##  semilogy(x,n6(ii, 1:N),'LineWidth', 0.2)
##  %endfor
##  title (sprintf ("%g", t(ii)));
##  legend ('e');
##  %axis ([min(x) max(x) 0 max(max(state0))]);
##  subplot(1, 2, 2)
##  plot (t, V, t(ii), V(ii), 'ro')
##  title ('V')
##  xlabel ('t [s]')
##  ylabel ('V [V]')
##  drawnow
##endfor

##figure()
##semilogy(x,n1(1, 1:N)*2.5e19,'LineWidth', 0.5,'ro')
##semilogy(x,n1(end, 1:N)*2.5e19,'LineWidth', 0.5,'rx')
##hold on
##semilogy(x,n2(1, 1:N)*2.5e19,'LineWidth', 0.5,'go')
##semilogy(x,n2(end, 1:N)*2.5e19,'LineWidth', 0.5,'gx')
##
##semilogy(x,n3(1, 1:N)*2.5e19,'LineWidth', 0.5,'bo')
##semilogy(x,n3(end, 1:N)*2.5e19,'LineWidth', 0.5,'bx')
##
##semilogy(x,n3(1, 1:N)*2.5e19,'LineWidth', 0.5,'bo')
##semilogy(x,n3(end, 1:N)*2.5e19,'LineWidth', 0.5,'bx')
##
##semilogy(x,n4(1, 1:N)*2.5e19,'LineWidth', 0.5,'co')
##semilogy(x,n4(end, 1:N)*2.5e19,'LineWidth', 0.5,'cx')
##
##semilogy(x,n5(1, 1:N)*2.5e19,'LineWidth', 0.5,'y0')
##semilogy(x,n5(end, 1:N)*2.5e19,'LineWidth', 0.5,'yx')
##
##
##semilogy(x,n6(1, 1:N)*2.5e19,'LineWidth', 0.5,'ko')
##semilogy(x,n6(end, 1:N)*2.5e19,'LineWidth', 0.5,'kx')
figure()
semilogy(x,n1(1, 1:N)*2.5e19,'LineWidth', 1,':')
hold on
semilogy(x,n2(1, 1:N)*2.5e19,'LineWidth', 1,'-.')
semilogy(x,n3(1, 1:N)*2.5e19,'LineWidth', 1, '--')
semilogy(x,n4(1, 1:N)*2.5e19,'LineWidth', 1, '-')
semilogy(x,n5(1, 1:N)*2.5e19,'LineWidth', 1, '.')
semilogy(x,n6(1, 1:N)*2.5e19,'LineWidth', 1)
legend (elements)
title('initialzation')
xlabel ('t [s]')
ylabel ('V [V]')
drawnow

figure()
semilogy(x,n1(end, 1:N)*2.5e19,'LineWidth', 1,':')
hold on
semilogy(x,n2(end, 1:N)*2.5e19,'LineWidth', 1, '-.')
semilogy(x,n3(end, 1:N)*2.5e19,'LineWidth', 1, '--')
semilogy(x,n4(end, 1:N)*2.5e19,'LineWidth', 1,'-')
semilogy(x,n5(end, 1:N)*2.5e19,'LineWidth', 1,'.')
semilogy(x,n6(end, 1:N)*2.5e19,'LineWidth', 1)
legend(elements)
title('final state')
xlabel ('t [s]')
ylabel ('V [V]')
drawnow
