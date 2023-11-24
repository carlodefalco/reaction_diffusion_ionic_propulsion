function dstate= compute_change_rates_Townsend(density, reactions, index,x, electron_flux)
  N=numel(x);
  N_species=numfields(index);
  state=reshape(density,N,N_species);
  dstate=zeros(size(density));
  state(1,:)=[];
  sum_Ri = zeros (size (state));

  for the_reaction = reactions(:)'


    Rfi =ones(N-1,1).* the_reaction.rate_coeffs(1);

    for [value, key] = the_reaction.reactants
      %xk=state(:, index.(key))./2.5e25;%% dasistemare
       xk=1;
      Rfi .*= xk.*abs(electron_flux);
    endfor

    Ri = Rfi;% - Rbi;

    for [value, key] = the_reaction.reactants
      sum_Ri (:,index.(key)) -= Ri.* value;

    endfor
    for [value, key] = the_reaction.products
      sum_Ri (:,index.(key)) += Ri.* value;

    endfor

    endfor
 dstate=[zeros(1,N_species);sum_Ri];

 endfunction
