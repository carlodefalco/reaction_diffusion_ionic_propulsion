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
state0(idx.("Ar"))   = 2.0e+17;
state0(idx.("e"))    = 2.0e+17;
state0(idx.("Ar+"))  = 2.0e+17;
state0(idx.("Ar2+")) = 2.0e+10;
state0(idx.("Ar*"))  = 2.0e+10;

%%NOTE THAT THE MOBILITY SHOULD BE COMPUTED NOT IMPOSED: this is for a try
mobility = zeros (numfields (idx), 1);
mobility(idx.("Ar"))   = 1e-3;
mobility(idx.("e"))    = 1e-3;
mobility(idx.("Ar+"))  = 1e-3;
mobility(idx.("Ar2+")) = 1e-3;
mobility(idx.("Ar*"))  = 1e-3;

valence = zeros (numfields (idx), 1);
valence(idx.("Ar"))   = 0;
valence(idx.("e"))    = -1;
valence(idx.("Ar+"))  = 1;
valence(idx.("Ar2+")) = +1;
valence(idx.("Ar*"))  = 0;

elements= fieldnames(idx);
N_species=numfields(idx);
%% Problem Data
L = 2e-3;
N = 81;
x = linspace (0, L, N)';
T0 = 0.99e-4;
T1 = 1.02e-4;

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

y0dot=zeros(N*length(elements),1);


%% build the system
fun=@(t,y, ydot) drift_diffusion_reaction_system (t, y, ydot, r, idx, x, N, V0, q, epsilon, Vth, mobility, valence);
options = odeset('RelTol', 10.0^(-7), 'AbsTol', 10.0^(-7), 'normcontrol', 'on', 'refine', 5,...
                'Jacobian', @(t, y, ydot) system_jacobian(t, y, ydot, r, idx, x, N, V0, q, epsilon, Vth, mobility, valence));

%integration in time
[t, y] = ode15i (fun, [T0, T1], y0, y0dot, options);

%% rearrenge y vector to plot results by species
for k= 1:length(elements)
    species(k).density = y(:, 1+ N*(k-1):k*N);
endfor

V = V0(t')(2, :);

for ii = 1 : numel (t)
  figure()
  subplot(1, 2, 1)
  for k= 1:length(elements)
  plot (x,species(k).density(ii, 1:N),'LineWidth', 2)
  hold on
  endfor
  title (sprintf ("%g", t(ii)));
  legend (elements);
  axis ([min(x) max(x) 0 max(max(state0(1)))]);
  subplot(1, 2, 2)
  plot (t, V, t(ii), V(ii), 'ro')
  title ('V')
  xlabel ('t [s]')
  ylabel ('V [V]')
  drawnow
endfor
