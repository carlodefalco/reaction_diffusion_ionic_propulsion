
function dstate = compute_change_rates3 (density, reactions, index,x)
  N=numel(x);
  N_species=numfields(index);
  state=reshape(density,N,N_species);
  dstate=zeros(size(density));

  sum_Ri = zeros (size (state));
  for the_reaction = reactions(:)'


    Rfi =ones(N,1).* the_reaction.rate_coeffs(1);

    for [value, key] = the_reaction.reactants
      Rfi .*= state (:,index.(key)).^ value;
    endfor

    Rbi = ones(N,1).* the_reaction.rate_coeffs(2);
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
