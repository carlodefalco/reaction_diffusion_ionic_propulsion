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

clc
clear

addpath (canonicalize_file_name ("../../../data"));
addpath (canonicalize_file_name ("./funzioni chimica manuali"));
addpath (canonicalize_file_name ("../"));
addpath (canonicalize_file_name ("../IPROP_ARGON/READ_REACTIONS"));
[r, idx] = read_reactions (file_in_loadpath ("balcon_et_al_argon_ionization3eV.json"));
##[r, idx] = read_reactions (file_in_loadpath ("tamburini10Td.json"));
elements=fieldnames(idx);
pretty_print_reactions (r);

##
x0 = zeros (numfields (idx), 1); %by defaults the photon initial density is zero.
x0(idx.("Ar"))   = 2.5e+25;
x0(idx.("e"))    = 1.01e19;
x0(idx.("Ar+"))  = 1.0e19;
x0(idx.("Ar2+")) = 1.0e17;
x0(idx.("Ar*"))  = 1.0e20;
x0(idx.("h_nu")) = 0;

T0   = 0;
Tend = 1e-3;
T_vec=[0 logspace(-10, log10 (Tend), 1000)];
s_e0=x0(idx.("e"));
s_Ar0=x0(idx.("Ar"));
s_Ar_excited0=x0(idx.("Ar*"));
s_Ar_plus0=x0(idx.("Ar+"));
s_Ar2_plus0=x0(idx.("Ar2+"));
s_hnu0=x0(idx.("h_nu"));
R1_f0 = r(1).rate_coeffs(1)*s_e0*s_Ar0;
R2_f0 = r(2).rate_coeffs(1)*s_e0*s_Ar0;
R3_f0 = r(3).rate_coeffs(1)*s_e0*s_Ar_excited0;
R4_f0 = r(4).rate_coeffs(1)*s_Ar_excited0^2;
R5_f0 = r(5).rate_coeffs(1)*s_Ar_plus0*s_Ar0^2;
R6_f0 = r(6).rate_coeffs(1)*s_e0*s_Ar2_plus0;
R7_f0 = r(7).rate_coeffs(1)*s_Ar_excited0;
R8_f0 = r(8).rate_coeffs(1)*s_e0*s_Ar0;
s_e0_dot = R1_f0 + R3_f0 + R4_f0 - R6_f0;
s_Ar0_dot = - R1_f0 - R2_f0 + R4_f0 -R5_f0 + R6_f0 + R7_f0;
s_Ar_plus0_dot = R1_f0 + R3_f0 + R4_f0 - R5_f0;
s_Ar2_plus0_dot= + R5_f0 - R6_f0;
s_Ar_excited0_dot= + R2_f0 - R3_f0- 2*R4_f0 + R6_f0 - R7_f0;
s_hnu0_dot= R7_f0;

stiff_par=([x0(idx.("e")) * s_e0_dot ;
x0(idx.("Ar")) *s_Ar0_dot;
x0(idx.("Ar+")) *s_Ar_plus0_dot ;
x0(idx.("Ar2+")) *s_Ar2_plus0_dot;
x0(idx.("Ar*")) *s_Ar_excited0_dot;
x0(idx.("h_nu")) *s_hnu0_dot]).^-1;
M=eye(6);

x0dot=[s_e0_dot, s_Ar0_dot,s_Ar_plus0_dot,s_Ar_excited0_dot,s_Ar2_plus0_dot,s_hnu0_dot];




fun_jac = @(x) (- compute_change_rates (x, r, idx));

options3=odeset('Jacobian',  @(t, x, xdot) implicit_change_rates_jacobian(t,x,xdot,r, idx));


%%Integration

eqs = @(t, x, xdot) compute_change_rates_implicit(x, xdot, r, idx);
[t, y] = ode15i ( eqs, T_vec, x0, x0dot, options3);


figure()
for [val, key] = idx
  warning ("off", "Octave:negative-data-log-axis")
  loglog (t(2:end), y(:, val)(2:end), 'Linewidth', 2)
  set (gca, 'fontsize', 24,
       'ytick', logspace(1, 20, 9),
       'xtick', logspace(-10, 5, 6))
  ylim ([1e1 1e26])
  title('ode15i')
  hold all
endfor
legend (fieldnames (idx){:}, 'location', 'eastoutside')

##figure()
##
##for [val, key] = idx
##  warning ("off", "Octave:negative-data-log-axis")
##  loglog (t3(2:end), y3(:, val)(2:end), 'Linewidth', 2)
##  set (gca, 'fontsize', 24,
##       'ytick', logspace(-12, 20, 9),
##       'xtick', logspace(-15, 5, 6))
##  ylim ([1e-12 1e26])
##   title('euler')
##  hold all
##endfor
##legend (fieldnames (idx){:}, 'location', 'eastoutside')
##hold off
