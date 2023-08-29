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
[r, idx] = read_reactions (file_in_loadpath ("robertson_autocatalysis.json"));

pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1);
x0(idx.("A"))  = 2.5e19;
x0(idx.("B"))  = 1e6;
x0(idx.("C"))  = 1e6;

T0   = 0;
Tend = 1e6;

o = odeset ('RelTol',1e-4, 'AbsTol', [1e-6 1e-10 1e-6], 'Jacobian', @(t, x)  compute_change_rates_jacobian (x, r, idx));
[t, y] = ode15s (@(t, x) compute_change_rates (x, r, idx),[0 4*logspace(-6,6)], x0, o);



figure
semilogx (t, y(:,1))
hold all
semilogx (t, 1e4*y(:,2))
semilogx (t, y(:,3))
legend ("A","B*1e4","C")

hold off
