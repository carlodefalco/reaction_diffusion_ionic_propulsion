
%% Solve :
%%   n_k' = div (mun Vth (grad n_k - grad phi/Vth n_k))
%%   0  = -div (epsilon grad (phi)) - q sum_k { z_k  n_k}
%%   n_k(0) = n_k(L)= 0
%%   phi(0) = 0 phi(L) = V0(t)
%%   this function builds the jacobian of the system of equation
function [Jac, DJac]= system_jacobian (t, y, ydot, reactions, index, x, N, V0, q, epsilon, Vth, mobility, valence)

  V0t = V0(t);
  N_species= numfields(index);
  rearrange_y=zeros(N, N_species);
  for k=1:N_species
    for ii= 1:N
        rearrange_y(ii, k)= y(ii+(k-1)*N);
    endfor
   endfor
  Jabian_reaction = sparse (N_species*N, N_species*N);

  for ii=1:N
    Jacobian_reaction(ii:N:ii+N*(N_species-1), ii:N:ii+N*(N_species-1)) = compute_change_rates_jacobian(rearrange_y(ii,:), reactions, index);
  endfor

  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  sum_k=zeros(N,1);
  for k=1:N_species
    sum_k+= rearrange_y(:,k).*valence(k);
  endfor
  b=q*M*sum_k;
  phi = zeros (N, 1);
  phi([1 N]) = V0t;

  phi(2:end-1) = A00(2:end-1, 2:end-1) \ (b(2:end-1)-A00(2:end-1, :) * phi);


  J = sparse (N_species*N, N_species*N);
  for k=1:N_species
    An = - bim1a_advection_diffusion (x, mobility(k)*Vth, 1, 1, -valence(k)*phi/Vth);
    An([1 end], :) = 0;
    An([1 end], [1 end]) = 1;

    J(1+(N*(k-1)): k*N, 1+(N*(k-1)):k*N) = M\An;
  endfor
  DJac = sparse(N*N_species, N*N_species);
  for k=1:N_species
    Mass(1+(N*(k-1)): k*N, 1+(N*(k-1)):k*N) = M;
  endfor
  DJac=sparse(Mass);
  Jac=(J+Jacobian_reaction);
endfunction
