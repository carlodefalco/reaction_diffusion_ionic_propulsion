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
addpath (canonicalize_file_name ("../"));
[r, idx] = read_reactions (file_in_loadpath ("balcon_et_al_argon_ionization.json"));

pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1); %by defaults the photon initial density is zero.

x0(idx.("Ar"))   = 2.5e+19/2.5e+19;
x0(idx.("e"))    = 1.0e+6/2.5e+19;
x0(idx.("Ar+"))  = 1.0e+6/2.5e+19;
x0(idx.("Ar2+")) = 1.0e+3/2.5e+19;
x0(idx.("Ar*"))  = 1.0e+10/2.5e+19;


x0dot=zeros(numfields (idx), 1);
T0   = 0;
Tend = 1.0e-5;
T=linspace(T0, Tend, 200);

eqs = @(t, x, xdot) compute_change_rates_implicit (x, xdot, r, idx);
options = odeset('RelTol', 10.0^(-7), 'AbsTol', 10.0^(-7), 'Jacobian', @(t, x, xdot) implicit_change_rates_jacobian(t, x, xdot, r, idx));
[t, y] = ode15i ( eqs, T, x0, x0dot, options);

y=y.*2.5e+19;
figure
for ii=1:numfields(idx)
loglog (T,y(:,ii), 'LineWidth', 2)

hold on
endfor

legend (fieldnames (idx){:}, 'location', 'eastoutside', 'FontSize', 18)
hold off
