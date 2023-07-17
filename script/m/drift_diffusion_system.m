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

species(1).id= 'e';
species(1).mobility=1e-3;
species(1).valence=-1;
species(1).n0=[];
species(1).density=[];

species(2).id= 'p';
species(2).mobility=1e-4;
species(2).valence=1;
species(2).n0=[];
species(2).density=[];

species(3).id= 'Ar+';
species(3).mobility=1e-5;
species(3).valence=1;
species(3).n0=[];
species(3).density=[];





L = 2e-3;
N = 81;
x = linspace (0, L, N).';
T0 = 0.999e-4;
T1 = 1.02e-4;
T_vec=linspace (T0, T1, 40);

V0 = @(t) 5000*[zeros(size(t)); ((t-1e-4)*1e6 .* ((t>=1e-4)&(t<=1.01e-4)) + 1.0 .* (t>1.01e-4))];
q  = 1.6e-19; %elementary charge
epsilon = 8.8e-12; %dielectric constant in vacuum
%% Diffusion coefficient according to Einstein–Smoluchowski relations is D_k= mu_k Vth
Vth = 26e-3;

N_species=size(species,2);
y0=[];


%for k= 1:N_species
%  species(k).n0= x .* (L - x)*k * 2e17 / (L/2)^2;
%  y0=[y0; species(k).n0];
%endfor
species(1).n0 = x .* (L - x) * 2e17 / (L/2)^2;
species(2).n0 = x .* (L - x) * 2e17 / (L/2)^2;
species(3).n0 = x .* (L - x) * 2e15 / (L/2)^2;
for k= 1:N_species
  y0=[y0; species(k).n0];
endfor
o = odeset ('Jacobian', @(t, y) jacobian (t, y, x, N, N_species, V0, q, epsilon, Vth, species),
	    'InitialStep', 1e-4-T0);

[t, y] = ode15s (@(t, y) odefun (t, y, x, N, N_species, V0, q, epsilon, Vth, species),
		 T_vec, y0,o);


for k= 1:N_species
    species(k).density = y(:, 1+ N*(k-1):k*N);
endfor

V = V0(t')(2, :);

for ii = 1 : numel (t)
  subplot(1, 2, 1)
  plot (x,species(1).density(ii, 1:N), x,species(2).density(ii, 1:N), x,species(3).density(ii, 1:N))
  title (sprintf ("%g", t(ii)));
  legend ('e', 'p', 'Ar');
  axis ([min(x) max(x) 0 max(max(species(1).n0))]);
  subplot(1, 2, 2)
  plot (t, V, t(ii), V(ii), 'ro')
  title ('V')
  xlabel ('t [s]')
  ylabel ('V [V]')
  drawnow
endfor
