%% Solve :
%%   n_k' = div (mun Vth (grad n_k - grad phi/Vth n_k))
%%   0  = -div (epsilon grad (phi)) - q sum_k { z_k  n_k}
%%   n_k(0) = n_k(L)= 0
%%   phi(0) = 0 phi(L) = V0(t)
%%   this function builds the jacobian of the system of equation
function J = jacobian (t, y, x, N, N_species, V0, q, epsilon, Vth, species)

  V0t = V0(t);

  for k=1:N_species
    species(k).density=y(1+(N*(k-1)):N*k);
  endfor

  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  a=0;
  for k=1:N_species
    a+= species(k).density*species(k).valence;
  endfor
  b=q*M*a;
  phi = zeros (N, 1);
  phi([1 N]) = V0t;

  phi(2:end-1) = A00(2:end-1, 2:end-1) \ (b(2:end-1)-A00(2:end-1, :) * phi);


  J = sparse (N_species*N, N_species*N);
  for k=1:N_species
    An=zeros(N,N);
    An = - bim1a_advection_diffusion (x, species(k).mobility*Vth, 1, 1, -species(k).valence*phi/Vth);
    An([1 end], :) = 0;
    An([1 end], [1 end]) = 1;

    J(1+(N*(k-1)): k*N, 1+(N*(k-1)):k*N) = M\An;
  endfor

endfunction
