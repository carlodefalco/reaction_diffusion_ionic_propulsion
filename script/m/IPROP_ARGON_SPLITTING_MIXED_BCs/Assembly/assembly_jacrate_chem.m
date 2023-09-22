%function [jacrate]=assembly_jacrate_chem(t, y, x, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk,ui, tau)
function [jacrate]=assembly_jacrate_chem(t, y, x,reactions, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk)
  N_species=numfields(index);
  M=kron(eye(N_species), Mk);
  %V0t=bc_phi(t);
  %y_bc=y;
%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================

  %p=y_bc(1:N-2);
  %n=y_bc(1+(N-2):2*(N-2));
  %dsource_p=diag((-n.^2.-ui^2)./((p+n).^2.*tau));
  %dsource_n=diag((-p.^2.-ui^2)./((p+n).^2.*tau));
  %DF=sparse([Mk(2:end-1, 2:end-1)*dsource_p, Mk(2:end-1, 2:end-1)*dsource_n; Mk(2:end-1, 2:end-1)*dsource_p, Mk(2:end-1, 2:end-1)*dsource_n]);

##  density=M*y_bc(1:(N-2)*N_species);
##  DF=compute_change_rates_jacobian3(density, reactions, index,x(2:end-1)); %jacobian matrix of the chemistry(N-2)* N_species x (N-2)* N_species
  density=M*y(1:N*N_species);
  DF=compute_change_rates_jacobian3(density, reactions, index,x); %jacobian matrix of the chemistry(N-2)* N_species x (N-2)* N_species

  %jac_chem=sparse((N_species)*(N-2),(N_species)*(N-2));
  jac_chem=sparse(DF);
  jacrate=+jac_chem;

endfunction
