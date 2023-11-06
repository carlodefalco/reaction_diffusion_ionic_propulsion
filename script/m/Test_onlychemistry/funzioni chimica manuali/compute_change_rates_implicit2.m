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


## rates = compute_change_rates (state, reactions, index)
##
## given the set of M reactions involving N species
##
## \sum_{k=1}^{N} r_{i,k} s_k <--> \sum_{k=1}^{N} p_{i,k} s_k    i=1..M
##
## with forward reaction coefficient cf_i and backward reaction coefficient cb_i
##
## the rate Rf_i of the i-th forward reaction is
##
## Rf_i = cf_i \prod_k s_k^{\left( r_{i,k} \right)}
##
## the rate Rb_i of the i-th forward reaction is
##
## Rb_i = cb_i \prod_k s_k^{\left( p_{i,k} \right)}
##
## the net rate of the i-th reaction is
##
## R_i = Rf_i - Rb_i
##
## the system of odes is
##
## \dot{s}_k = \sum_{i=1}^M \left( - Rf_i r_{i,k} + Rf_i p_{i,k}
##             - Rb_i p_{i,k} - Rb_i r_{i,k} \right) =
##             - \sum_{i=1}^M \left( R_i (r_{i,k} - p_{i,k}) \right)





function Res = compute_change_rates_implicit2 (state, dstate, reactions, index)

  s_A=state(index.("A"));
  s_B=state(index.("B"));
  s_C=state(index.("C"));

  s_A_dot = dstate(index.("A"));
  s_B_dot = dstate(index.("B"));
  s_C_dot = dstate(index.("C"));

  R1_f = reactions(1).rate_coeffs(1)*s_A;
  R2_f = reactions(2).rate_coeffs(1)*s_B^2;
  R3_f = reactions(3).rate_coeffs(1)*s_B*s_C;

  sum_Ri_A = -R1_f + R3_f;
  sum_Ri_B = R1_f - R2_f - R3_f;
  sum_Ri_C = R2_f;

  Res=[ s_A_dot - sum_Ri_A;
        s_B_dot - sum_Ri_B;
        s_C_dot - sum_Ri_C
       ];




endfunction
