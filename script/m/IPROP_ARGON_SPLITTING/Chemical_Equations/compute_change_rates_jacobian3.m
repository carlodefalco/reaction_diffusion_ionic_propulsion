function jac = compute_change_rates_jacobian3( density, reactions, index,x)

  N_species=numfields(index);
  N=numel(x);
  state=reshape(density,N,N_species);
  jac = sparse(N*N_species, N*N_species);

  for the_reaction = reactions(:)'

  d_Rfi_d_sj=zeros(N,N_species);

   for [value, key] = the_reaction.reactants
      d_Rfi_d_sj(:,index.(key)) = ones(N,1).*the_reaction.rate_coeffs(1) .* state(:,index.(key)).^(value-1)* value;% derivative with respect to the index.(key) concentration
      for [value2, key2] = the_reaction.reactants
        if index.(key2)!=index.(key);
          d_Rfi_d_sj(:,index.(key)).*=state(:,index.(key2)).^(value2); % multiply by the other concentrations wrt no derivative is taken
        endif
      endfor
    endfor

   d_Rbi_d_sj=zeros(N,N_species);
##   for [value, key] = the_reaction.products
##      d_Rbi_d_sj(index.(key)) = the_reaction.rate_coeffs(2) * state(index.(key)) ^ (value-1)* value;
##      for [value2, key2] = the_reaction.products
##        if index.(key2)!=index.(key);
##          d_Rbi_d_sj(index.(key))*=state(index.(key2))^(value2);
##        else
##          d_Rbi_d_sj(index.(key))= d_Rbi_d_sj(index.(key));
##        endif
##      endfor
##    endfor

    d_Ri_d_sj = d_Rfi_d_sj- d_Rbi_d_sj;
    dstate_jac=sparse(N, N_species*N);
    for k=1:N_species
      dstate_jac(:, 1+N*(k-1):N*k)=diag(d_Ri_d_sj(:,k));
    endfor

    for [value, key] = the_reaction.reactants
      jac ( (index.(key)-1)*N+1:index.(key)*N,:) -= dstate_jac.* value;
    endfor
    for [value, key] = the_reaction.products
      jac ((index.(key)-1)*N+1:index.(key)*N,:) += dstate_jac.* value;
    endfor
  endfor




endfunction
