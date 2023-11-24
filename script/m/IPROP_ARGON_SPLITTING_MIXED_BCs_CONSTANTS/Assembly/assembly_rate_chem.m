function [rate]=assembly_rate_chem( y, x, reactions, index, valence, mobility,N, q, epsilon, Vth, P, Mk,electron_flux)
 N_species=numfields(index);
 %V0t=bc_phi;
%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================

  chemistry=compute_change_rates_Townsend(y(1:N*N_species), reactions, index,x, electron_flux);

## chemistry= compute_change_rates3 (y(1:N*N_species), reactions, index, x); %chemical source term ((N-2)x N_species)
## F= zeros(N*N_species,1); %compute source term wihtout boundary conditions
## for k=1:N_species
##  F(1+N*(k-1):k*N)  = Mk* chemistry (:,k);
## endfor
f1=sparse(N*N_species,1);
 for k=1:N_species
    F=bim1a_rhs(x,chemistry(:,k),1);
    f1(1+(k-1)*N:k*N)=F;
   endfor
 rate =f1;

endfunction
