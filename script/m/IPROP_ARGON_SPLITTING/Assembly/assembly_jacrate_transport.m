%function [jacrate]=assembly_jacrate_transport(t, y, x, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk,ui, tau)
function [jacrate]=assembly_jacrate_transport(t,y, x,reactions, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk)
  N_species=numfields(index);
  M=kron(eye(N_species), Mk(2:end-1, 2:end-1));
  V0t=bc_phi(t);
  y_bc=y;

  A=sparse(N_species*N, N_species*N ); %global transport matrix (N_species*N x N_species*N )
  A_bc=sparse(N_species*(N-2),N_species*(N-2));%remove the boundary raws and colums N_species*(N-2) x N_species*(N-2)
  phi=[V0t(1); y_bc((N-2)*(N_species)+1:end); V0t(2) ];
  for k=1:N_species
    A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, -valence(k)*phi/Vth(k));
    A_bc(1+(N-2)*(k-1):(N-2)*k,1+(N-2)*(k-1):(N-2)*k) = A(2+N*(k-1):N*k-1,2+N*(k-1):N*k-1);
  endfor

%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================


  L=sparse([(A_bc), zeros(N_species*(N-2), N-2);
              kron(valence, -q.*Mk(2:end-1, 2:end-1))', P(2:end-1, 2:end-1)]);%left-hand side matrix

  jacrate=-L;

endfunction
