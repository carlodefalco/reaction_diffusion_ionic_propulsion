%function [rate]=assembly_rate_chem(t, y, x, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk,ui, tau)
function [rate]=assembly_rate_chem(t, y, x, reactions, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk)
 N_species=numfields(index);
 V0t=bc_phi(t);
 y_bc=y;
 coeff=1e22;
 D= coeff.*(x<x(end)/2);
 C= coeff.*(x>=x(end)/2);
%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================

 %p=y_bc(1:N-2);
 %n=y_bc(1+(N-2):2*(N-2));
 %source=(ui^2.-p.*n)./((p+n) .*tau);
 %chemistry=[source, source];
 chemistry= compute_change_rates3 (y_bc(1:(N-2)*N_species), reactions, index, x(2:end-1)); %chemical source term ((N-2)x N_species)
 F= zeros((N-2)*N_species,1); %compute source term wihtout boundary conditions
 for k=1:N_species
  F(1+(N-2)*(k-1):k*(N-2))  = Mk(2:end-1, 2:end-1)* chemistry (:,k);
 endfor

## f2=sparse([F; zeros(size(P(2:end-1,1)))]);
f2=sparse(F);
 rate = +f2;

 endfunction
