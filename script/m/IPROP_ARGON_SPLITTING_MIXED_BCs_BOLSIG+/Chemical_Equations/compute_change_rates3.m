function dstate = compute_change_rates3 ( variables, reactions, index,x)
  N=numel(x);
  N_species=numfields(index);
  density=variables(1:N*N_species);
  phi=variables(N_species*N+1:end);
  state=reshape(density,N,N_species);
  dstate=zeros(size(density));


  r=reactions;

  sum_Ri = zeros (size (state));
  for the_reaction = r(:)'

    Rfi = the_reaction.rate_coeffs(:,1);%ones(N,1).*

    for [value, key] = the_reaction.reactants
      Rfi .*= state (:,index.(key)).^ value;
    endfor

    Rbi = the_reaction.rate_coeffs(:,2);
    for [value, key] = the_reaction.products
      Rbi .*= state (:,index.(key)).^ value;
    endfor

    Ri = Rfi - Rbi;

    for [value, key] = the_reaction.reactants
      sum_Ri (:,index.(key)) -= Ri.* value;

    endfor
    for [value, key] = the_reaction.products
      sum_Ri (:,index.(key)) += Ri.* value;

    endfor

    endfor
 dstate=sum_Ri;




endfunction
