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



function dstate = compute_change_rates (state, reactions, index)

  dstate = zeros (size (state));
  for the_reaction = reactions(:)'

    Rfi = the_reaction.rate_coeffs(1);
    for [value, key] = the_reaction.reactants
      Rfi *= state (index.(key)) ^ value;
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
