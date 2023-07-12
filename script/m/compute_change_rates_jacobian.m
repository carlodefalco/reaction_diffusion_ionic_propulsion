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
##
## the jacobian of \dot{s}_k is given by
##
## \dfrac{\partial \dot{s}_k}{\partial s_j} = (r_{i,k} - p_{i,k}) \dfrac{\partial R_i}{\partial s_j}


function dstate_jac = compute_change_rates_jacobian (state, reactions, index)

  N = numel (state);
  dstate_jac = sparse (N, N);
  
  for the_reaction = reactions(:)'

    d_Rfi_d_sj = zeros (1, N);
    for [value, key] = the_reaction.reactants
      d_Rfi_d_sj(index.(key)) = the_reaction.rate_coeffs(1) * state (index.(key)) ^ (value-1) * value;
    endfor

    Rbi = the_reaction.rate_coeffs(2);
    for [value, key] = the_reaction.products
      Rbi *= state (index.(key)) ^ value;
    endfor

    Ri = Rfi - Rbi;

    for [value, key] = the_reaction.reactants
      dstate (index.(key)) -= Ri * value;
    endfor
    for [value, key] = the_reaction.products
      dstate (index.(key)) += Ri * value;
    endfor
    
  endfor

endfunction
