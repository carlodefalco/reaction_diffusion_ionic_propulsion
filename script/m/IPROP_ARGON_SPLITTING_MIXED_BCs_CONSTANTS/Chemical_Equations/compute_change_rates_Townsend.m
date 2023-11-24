function dstate= compute_change_rates_Townsend(density, reactions, index,x, electron_flux)
  N=numel(x);
  q=1.6e-19;
  N_species=numfields(index);
  state=reshape(density,N,N_species);
  dstate=zeros(size(density));
  state(1,:)=[];
  sum_Ri = zeros (size (state));

  for the_reaction = reactions(:)'

    target_fraction=1;
    Rfi =57.*target_fraction.*(abs(electron_flux)./q);
    Ri = Rfi;% - Rbi;

    for [value, key] = the_reaction.reactants
      sum_Ri (:,index.(key)) -= Ri.* value;

    endfor
    for [value, key] = the_reaction.products
      sum_Ri (:,index.(key)) += Ri.* value;

    endfor

    endfor
 dstate=sum_Ri;

 endfunction
