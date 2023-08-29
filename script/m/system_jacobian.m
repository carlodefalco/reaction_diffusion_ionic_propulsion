
%% Solve :
%%   n_k' = div (mun Vth (grad n_k - grad phi/Vth n_k))
%%   0  = -div (epsilon grad (phi)) - q sum_k { z_k  n_k}
%%   n_k(0) = n_k(L)= 0
%%   phi(0) = 0 phi(L) = V0(t)
%%   this function builds the jacobian of the system of equation
function Jac = system_jacobian (t, y, reactions, index, x, N, V0, q, epsilon, Vth, mobility, valence)

  V0t = V0(t);
  N_species= numfields(index);
  rearrange_y=zeros(N, N_species);
  for k=1:numel(x)
    for ii=1:N_species
      rearrange_y(k, ii)= y(k+(ii-1)*N);
    endfor
  endfor
  Jabian_reaction = sparse (N_species*N, N_species*N);
  for k=1:numel(x)
   jacobian_local=compute_change_rates_jacobian(rearrange_y(k,:),reactions,index);
   for ii=1:N_species
  Jacobian_reaction(k:N:k+N*(N_species-1),k:N:k+N*(N_species-1))=jacobian_local;
   endfor
  endfor

  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  a=zeros(N,1);
  for k=1:N_species
    a+= rearrange_y(:,k).*valence(k);
  endfor
  b=q*M*a;
  phi = zeros (N, 1);
  phi([1 N]) = V0t;

  phi(2:end-1) = A00(2:end-1, 2:end-1) \ (b(2:end-1)-A00(2:end-1, :) * phi);


  J = sparse (N_species*N, N_species*N);
  for k=1:N_species
    An=zeros(N,N);
    An = - bim1a_advection_diffusion (x, mobility(k)*Vth, 1, 1, -valence(k)*phi/Vth);
    An([1 end], :) = 0;
    An([1 end], [1 end]) = 1;

    J(1+(N*(k-1)): k*N, 1+(N*(k-1)):k*N) = M\An;
  endfor
  Jac=+J;
endfunction
