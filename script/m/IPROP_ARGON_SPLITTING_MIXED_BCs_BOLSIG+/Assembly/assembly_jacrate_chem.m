function [jacrate]=assembly_jacrate_chem(y, x,reactions, index, valence, mobility,bc, bc_phi,N,N_species, q, epsilon, Vth, P, Mk)
  N_species=numfields(index);
  M=kron(eye(N_species), Mk);
  %V0t=bc_phi(t);
%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================

  DF=compute_change_rates_jacobian3( y(1:N*N_species), reactions, index,x); %jacobian matrix of the chemistry(N-2)* N_species x (N-2)* N_species

  jac_chem=sparse(DF);
  jacrate=+jac_chem;

endfunction
