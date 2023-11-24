function [rate]=assembly_rate_transport(y, x, reactions, index, valence, mobility, q, epsilon, Vth, P, Mk)
 N_species=numfields(index);
 M=kron(eye(N_species),Mk);
 N=numel(x);
## for [value, key]=index
##  if (index.(key)== "h_nu")
##      N_species-=1;
##  endif
## endfor
  %===============================================================================
  % RATE
  %===============================================================================

  A=sparse(N_species*N, N_species*N ); %global transport matrix (N_species*N x N_species*N )
  phi=y(N_species*N+1:end);
  for k=1:N_species
        A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, -valence(k).*phi./Vth(k));%-valence(k)*phi/Vth(k)
  endfor

  K=kron(valence, -q.*Mk)';
  L=sparse([A, zeros(N_species*N, N);
            K , P]);%left-hand side matrix
  rate = -L*y;

endfunction
