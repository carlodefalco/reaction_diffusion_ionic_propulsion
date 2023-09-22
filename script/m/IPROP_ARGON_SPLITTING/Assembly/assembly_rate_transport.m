%function [rate]=assembly_rate_transport(t, y, x, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk,ui, tau)
function [rate]=assembly_rate_transport(t, y, x, reactions, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk)
 N_species=numfields(index);
 V0t=bc_phi(t);
 y_bc=y;
 coeff=1e22;
 D= coeff.*(x<x(end)/2);
 C= coeff.*(x>=x(end)/2);
%===============================================================================
% RATE
%===============================================================================

 A=sparse(N_species*N, N_species*N ); %global transport matrix (N_species*N x N_species*N )
 A_bc=sparse(N_species*(N-2),N_species*(N-2));%remove the boundary raws and colums N_species*(N-2) x N_species*(N-2)
 phi=[V0t(1); y_bc((N-2)*(N_species)+1:end); V0t(2) ];
 for k=1:N_species
   A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, -valence(k)*phi/Vth(k));
   A_bc(1+(N-2)*(k-1):(N-2)*k,1+(N-2)*(k-1):(N-2)*k) = A(2+N*(k-1):N*k-1,2+N*(k-1):N*k-1);
 endfor

 boundary_RHS=zeros((N-2)*N_species,1); %boundary conditions contribution to the righ hand side
 for k=1:N_species
  boundary_RHS(1+(N-2)*(k-1):k*(N-2))= - A(2+N*(k-1):N*k-1, N*(k-1)+1).*bc(k,1) - A(2+N*(k-1):N*k-1,N*k).*bc(k,2);
 endfor
 L=sparse([A_bc, zeros(N_species*(N-2), N-2);
           kron(valence, -q.*Mk(2:end-1, 2:end-1))', P(2:end-1, 2:end-1)]);%left-hand side matrix


 %f1=sparse([boundary_RHS; -P(2:end-1,1).*V0t(1)-P(2:end-1,end).*V0t(2)+(q.*Mk*(D-C))(2:end-1)]); %right hand side matrix -DF
 f1=sparse([boundary_RHS; -P(2:end-1,1).*V0t(1)-P(2:end-1,end).*V0t(2)]); %right hand side matrix -DF

 rate = -L*y_bc+f1;


 endfunction
