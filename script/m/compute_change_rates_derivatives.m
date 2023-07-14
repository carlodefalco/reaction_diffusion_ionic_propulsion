function  d_Rfi_d_sj = compute_change_rates_derivatives(state, reactions, index)

  N = numel (state);
  d_Rfi_d_sj=sparse(N,N);
  for the_reaction = reactions(:)'
    for [value, key] = the_reaction.reactants
      d_Rfi_d_sj(index.(key)) = the_reaction.rate_coeffs(1) * state(index.(key)) ^ (value-1)* value;
      for [value2, key2] = the_reaction.reactants
        if key2!=key;
          d_Rfi_d_sj(index.(key))*=state(index.(key2))^(value2);
        else
          d_Rfi_d_sj(index.(key))= d_Rfi_d_sj(index.(key));
        endif
      endfor
    endfor
   endfor
 endfunction
