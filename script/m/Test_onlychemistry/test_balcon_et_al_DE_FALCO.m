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


addpath (canonicalize_file_name ("../../../data"));
addpath (canonicalize_file_name ("./funzioni chimica manuali"));
addpath (canonicalize_file_name ("../"));
addpath (canonicalize_file_name ("../IPROP_ARGON/READ_REACTIONS"));
[r, idx] = read_reactions (file_in_loadpath ("balcon_et_al_argon_ionization.json"));

pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1);
x0(idx.("Ar"))   = 2.5e+25;
x0(idx.("e"))    = 1.01e+18;
x0(idx.("Ar+"))  = 1e18;
x0(idx.("Ar2+")) = 1e16;
x0(idx.("Ar*"))  = 1e18;
x0(idx.("h_nu")) = 0;

T0   = 0;
Tend = 1;

M = eye (6);
xdot0 = M \ compute_change_rates (x0, r, idx);
fun = @(t, x, xdot) (M*xdot - compute_change_rates (x, r, idx));
function [jx, jxdot] = jacfun (t, x, xdot, M, r, idx)
  persistent t0 = 0
  if t > t0
    disp(t)
    t0 = t;
  endif
  jx = - compute_change_rates_jacobian (x, r, idx);
  jxdot = M;
endfunction

o = odeset ('AbsTol', 1e-9, 'RelTol', 1e-9, 'Jacobian', @(t, x, xdot) jacfun (t, x, xdot, M, r, idx));
[t, x] = ode15i (fun, [0 logspace(-10, log10 (Tend), 1000)], x0, xdot0);% o);
T_vec=[0 logspace(-10, log10 (Tend), 1000)];
dt=(T_vec(2:end)-T_vec(1:end-1))';
dt=dt*ones(1,6);
deltax0=min(abs((x(2:end,:)-x(1:end-1,:))./dt));

delta=(x(2:end,:)-x(1:end-1,:))./dt;

figure
for [val, key] = idx
  warning ("off", "Octave:negative-data-log-axis")
  loglog (t(2:end), x(:, val)(2:end), 'Linewidth', 2)
  set (gca, 'fontsize', 24,
       'ytick', logspace(-12, 20, 9),
       'xtick', logspace(-15, 5, 6))
  ylim ([1e-12 1e26])
  hold all
endfor

legend (fieldnames (idx){:}, 'location', 'eastoutside')
hold off
