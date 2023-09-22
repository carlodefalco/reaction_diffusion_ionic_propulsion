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
%%   n' - div (mun Vth (grad n - grad phi/Vth n))=0;
%%   p' - div (mun Vth (grad p + grad phi/Vth p))=0;
%%    -div (epsilon grad (phi))=  q * (p - n)
%%   n(0) = n(L) = p(0) = p(L) = 0
%%   phi(0) = 0 phi(L) = V0(t)
clc
clear
pkg load fpl bim msh

L = 2e-3;
N = 81;
x = linspace (0, L, N).';
T0 = 0.999e-4;
T1 = 1.02e-4;

function dy = odefun (t, y, x, N, V0, q, epsilon, Vth, mun, mup)

  V0t = V0(t);

  n   = y(1:N-2);
  p   = y((N-2)+1:2*(N-2));

  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  b   = q*(M(2:end-1,2:end-1)*(p - n));

  phi = A00(2:end-1, 2:end-1) \ (b-A00(2:end-1, end) * V0t(2));
  phi=[ V0t(1); phi; V0t(2)];

  An =  bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);

  dn = An(2:end-1, 2:end-1) * n;
  dn =- M(2:end-1, 2:end-1) \ dn;

  Ap =  bim1a_advection_diffusion (x, mup*Vth, 1, 1, -phi/Vth);

  dp = Ap(2:end-1, 2:end-1) * p;
  dp = -M(2:end-1, 2:end-1) \ dp;

  dy = [dn; dp];

endfunction

function J = jacobian (t, y, x, N, V0, q, epsilon, Vth, mun, mup)
  V0t = V0(t);

  n   = y(1:N-2);
  p   = y((N-2)+1:2*(N-2));

  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  b   = q*(M(2:end-1,2:end-1)*(p - n));

  phi = A00(2:end-1, 2:end-1) \ (b-A00(2:end-1, end) * V0t(2));
  phi=[ V0t(1); phi; V0t(2)];


  J = sparse (2*(N-2), 2*(N-2));

  An =  bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);
  Ap =  bim1a_advection_diffusion (x, mup*Vth, 1, 1, -phi/Vth);



  J(1:N-2, 1:N-2) =- M(2:end-1, 2:end-1)\An(2:end-1, 2:end-1);
  J((N-2)+1:2*(N-2), (N-2)+1:2*(N-2)) =- M(2:end-1, 2:end-1)\Ap(2:end-1, 2:end-1);

endfunction

V0 = @(t) 5000*[zeros(size(t)); ((t-1e-4)*1e6 .* ((t>=1e-4)&(t<=1.01e-4)) + 1.0 .* (t>1.01e-4))];
q  = 1.6e-19;
epsilon = 8.8e-12;
mun = 1e-3;
mup = 1e-3;
Vth = 26e-3;

n0 = x .* (L - x) * 2e17 / (L/2)^2;
p0 = x .* (L - x) * 2e17 / (L/2)^2;

y0 = [n0(2:end-1); p0(2:end-1)];
o = odeset ('Jacobian', @(t, y) jacobian (t, y, x, N, V0, q, epsilon, Vth, mun, mup),
	    'InitialStep', 1e-4-T0);
[t, y] = ode15s (@(t, y) odefun (t, y, x, N, V0, q, epsilon, Vth, mun, mup),
		 linspace (T0, T1, 100), y0, o);

n = [zeros(size(t)),y(:, 1:N-2),zeros(size(t))];
p = [zeros(size(t)),y(:, (N-2)+(1:N-2)),zeros(size(t))];
V = V0(t')(2, :);
for ii = 1 : numel (t)
  subplot(1, 2, 1)
  plot (x,n(ii, :), x, p(ii, :))
  title (sprintf ("%g", t(ii)));
  legend ('n', 'p');
  axis ([min(x) max(x) 0 max(max(n0))]);
  subplot(1, 2, 2)
  plot (t, V, t(ii), V(ii), 'ro')
  title ('V')
  xlabel ('t [s]')
  ylabel ('V [V]')
  drawnow
endfor
