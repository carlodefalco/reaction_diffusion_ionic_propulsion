function [jacrate]=assembly_jacrate_chem( y, x,reactions, index, valence, mobility,N, q, epsilon, Vth, P, Mk,current)
 N_species=numfields(index);
 M=kron(eye(N_species), Mk);
%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================

## DF=compute_change_rates_jacobian3(density, reactions, index,x); %jacobian matrix of the chemistry(N-2)*
 DF=compute_change_rates_jacobian3(y(1:N*N_species), reactions, index,x); %jacobian matrix of the chemistry(N-2)* N_species x (N-2)* N_species
 %DF((index.("Ar+")-1)*N+1,(1:N:N_species*N))=0;
 %DF((index.("Ar2+")-1)*N+1,(1:N:N_species*N))=0;

 jac_chem=sparse(DF);
 jacrate=+jac_chem;

endfunction
