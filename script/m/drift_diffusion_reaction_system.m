%% Solve :
%%   n_k' =- div (mun Vth (grad n_k + z_k grad phi/Vth n_k)) + ndot
%%   0  = -div (epsilon grad (phi)) - q sum_k { z_k  n_k}
%%   n_k(0) = n_k(L)= 0
%%   phi(0) = 0 phi(L) = V0(t)
%%   this function builds the system of ODE in implicit form

function eqs = drift_diffusion_reaction_system (t, y, ydot, r, index, x, N, V0, q, epsilon, Vth, mobility, valence)

    V0t = V0(t);
    N_species= numfields(index);

  %% The following lines are to handle the function compute_change_rates. This function takes as impute the N_species x1 state vector and
  %% compute the variation of the state due to chemical reaction in a point. To give to this function the correct input the y vector (that is N*N_species x 1)
  %% must be rewritten as a matrix Nx Nspecies. In this way it is possible to pass to the compute_change_rates function each line of this matrix, representing the local unkowns.
    rearrange_y=zeros(N, N_species);
    source_term=zeros(N, N_species);
    for k= 1:N
      for ii=1:N_species
        rearrange_y(k, ii)= y(k+(ii-1)*N);
      endfor
    endfor

  %% the local species production must be rewritten in a form compatible with the rearrange_y, namely we need a source_term of the same dimension.
  %% source_term is a matrix NxN_spcies containging the production chemical term for each species by columns in each N grid point
    for k= 1:N
     production_local=compute_change_rates(rearrange_y(k,:),r,index);
     for ii=1:N_species
      source_term(k, ii)= production_local(ii);
     endfor
    endfor

  %%(index.(elements{k})
    A00 = bim1a_laplacian (x, epsilon, 1);
    M   = bim1a_reaction (x, 1, 1);
    a=zeros(N,1);
    for k=1:N_species
      a+= rearrange_y(:,k).* valence(k);
    endfor
    b=q*M*a;

    phi = zeros (N, 1);
    phi([1 N]) = V0t;
    phi(2:end-1) = A00(2:end-1, 2:end-1) \ (b(2:end-1)-A00(2:end-1, :) * phi);


    dy=[];
    for k=1:N_species
      dn = zeros (N, 1);
      An = - bim1a_advection_diffusion (x, mobility(k)*Vth, 1, 1, -valence(k)*phi/Vth);
      dn(2:end-1) = An(2:end-1, :) * rearrange_y(:, k)+ M(2:end-1,:)*source_term(:, k);

      dn = M \ dn;
      dy = [dy; dn];
    endfor

    eqs=ydot-dy;
endfunction
