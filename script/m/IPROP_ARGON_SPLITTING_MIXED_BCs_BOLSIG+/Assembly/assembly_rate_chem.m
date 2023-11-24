function [rate]=assembly_rate_chem( y, x, reactions, index, valence, mobility,bc, bc_phi,N,N_species, q, epsilon, Vth, P, Mk)
 N_species=numfields(index);
 %V0t=bc_phi(t);
%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================

 chemistry= compute_change_rates3 (y(1:N*N_species), reactions, index, x); %chemical source term (Nx N_species)
## F= zeros(N*N_species,1); %compute source term wihtout boundary conditions
## for k=1:N_species
##  F(1+N*(k-1):k*N)  = Mk* chemistry (:,k);
## endfor
 F=reshape(chemistry, N_species*N, 1);
 f2=sparse(F);
 rate = +f2;
endfunction
