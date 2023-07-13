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


addpath (canonicalize_file_name ("../../data"));
[r, idx] = read_reactions (file_in_loadpath ("balcon_et_al_argon_ionization.json"));

pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1);
x0(idx.("Ar"))   = 2.5e+19;
x0(idx.("e"))    = 1.0e+6;
x0(idx.("Ar+"))  = 1.0e+6;
x0(idx.("Ar2+")) = 1.0e+3;
x0(idx.("Ar*"))  = 1.0e+10;

T0   = 0;
Tend = 1.0e-7;


o = odeset ('Jacobian', @(t, x)  compute_change_rates_jacobian(x, r, idx),'InitialStep', T0+1e-8);
[t, y] = ode15s (@(t, x) compute_change_rates (x, r, idx), [T0 Tend], x0, o);



figure
for [val, key] = idx
  semilogy (t(2:end), x(:, val)(2:end))
  hold all
endfor
legend (fieldnames (idx){:}, 'location', 'eastoutside')
hold off
