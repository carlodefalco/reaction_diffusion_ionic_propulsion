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
%%   n' = div (mun Vth (grad n - grad phi/Vth n))
%%   p' = div (mun Vth (grad p + grad phi/Vth p))
%%   0  = -div (epsilon grad (phi)) - q * (p - n)
%%   n(0) = n(L) = p(0) = p(L) = 0
%%   phi(0) = 0 phi(L) = V0(t)

pkg load fpl bim msh

L = 2e-3;
N = 31;
x = linspace (0, L, N).';
T = 1e-3;

function dy = odefun (t, y, x, N, V0, q, epsilon, Vth, mun, mup)

  V0t = V0(t);

  n   = y(1:N);
  p   = y(N+1:2*N);
  
  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  b   = q*(M*(p - n));
  phi = zeros (N, 1);
  phi([1 N]) = V0t;

  phi(2:end-1) = A00(2:end-1, 2:end-1) \ (-A00(2:end-1, :) * phi);
  
  An = - bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);
  dn = zeros (N, 1);
  dn(2:end-1) = An(2:end-1, :) * n;
  dn = M \ dn;

  Ap = - bim1a_advection_diffusion (x, mup*Vth, 1, 1, -phi/Vth);
  dp = zeros (N, 1);
  dp(2:end-1) = Ap(2:end-1, :) * p;
  dp = M \ dp;
  
  dy = [dn; dp];
  
endfunction 

function J = jacobian (t, y, x, N, V0, q, epsilon, Vth, mun, mup)

  V0t = V0(t);

  n   = y(1:N);
  p   = y(N+1:2*N);

  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  b   = q*(M*(p - n));
  phi = zeros (N, 1);
  phi([1 N]) = V0t;

  phi(2:end-1) = A00(2:end-1, 2:end-1) \ (-A00(2:end-1, :) * phi);
  
  
  J = sparse (2*N, 2*N);

  An = - bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);
  An([1 end], :) = 0;
  An([1 end], [1 end]) = 1;
  
  Ap = - bim1a_advection_diffusion (x, mup*Vth, 1, 1, -phi/Vth);
  Ap([1 end], :) = 0;
  Ap([1 end], [1 end]) = 1;

  
  J(1:N, 1:N) = M\An;
  J(N+1:2*N, N+1:2*N) = M\Ap;
  
endfunction

V0 = @(t) 1000*[0; ((t-1e-4)*1e4 .* ((t>=1e-4)&(t<=2e-4)) + 1.0 .* (t>2e-4))];
q  = 1.6e-19;
epsilon = 8.8e-12;
mun = 1e-3;
mup = 1e-3;
Vth = 26e-3;

n0 = x .* (L - x) * 1e6 / (L/2)^2;
p0 = x .* (L - x) * 1e6 / (L/2)^2;

y0 = [n0; p0];
o = odeset ('Jacobian', @(t, y) jacobian (t, y, x, N, V0, q, epsilon, Vth, mun, mup),
	    'InitialStep', 1e-7);
[t, y] = ode15s (@(t, y) odefun (t, y, x, N, V0, q, epsilon, Vth, mun, mup), [0 T], y0, o);

n = y(:, 1:N);
p = y(:, N+(1:N));
for ii = 1 : numel (t)
  plot (x,n(ii, :), x, p(ii, :))
  title (sprintf ("%g", t(ii)));
  legend ('n', 'p');
  axis ([min(x) max(x) 0 max(max(n0))]);
  drawnow
endfor
