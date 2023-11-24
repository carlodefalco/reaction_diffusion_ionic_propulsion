%function [jacrate]=assembly_jacrate(t, y, x, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk,ui, tau)
function [jacrate]=assembly_jacrate( y, x,reactions, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk, k_b, me,T_e)
 N_species=numfields(index);
 M=kron(eye(N_species), Mk);
% V0t=bc_phi(t);
 %y_bc=y;




   A=sparse(N_species*N, N_species*N ); %global transport matrix (N_species*N x N_species*N )
   %A_bc=sparse(N_species*(N-2),N_species*(N-2));%remove the boundary raws and colums N_species*(N-2) x N_species*(N-2)
   %phi=[V0t(1); y_bc((N-2)*(N_species)+1:end); V0t(2) ];
    phi=y(N_species*N+1:end);
   for k=1:N_species
     A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, -valence(k)*phi/Vth(k));
     %A_bc(1+(N-2)*(k-1):(N-2)*k,1+(N-2)*(k-1):(N-2)*k) = A(2+N*(k-1):N*k-1,2+N*(k-1):N*k-1);
   endfor
  A_bc=A;
##  A_bc((index.("Ar+")-1)*N+1 ,:)=0;
##  A_bc((index.("Ar+")-1)*N+1,(index.("Ar+")-1)*N+1)=1;
##  A_bc((index.("Ar2+")-1)*N+1,:)=0;
##  A_bc((index.("Ar2+")-1)*N+1,(index.("Ar2+")-1)*N+1)=1;
  A_bc(index.("e")*N,index.("e")*N )=A_bc(index.("e")*N,index.("e")*N)+ 1*(k_b*T_e/(2*pi*me))^0.5;;
##  P([1 end],:)=0;
##  P(1,1)=1;
##  P(end, end)=1;
%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================


   L=sparse([(A_bc), zeros(N_species*(N), N);
              kron(valence, -q.*Mk)', P]);%left-hand side matrix

   %p=y_bc(1:N-2);
   %n=y_bc(1+(N-2):2*(N-2));
   %dsource_p=diag((-n.^2.-ui^2)./((p+n).^2.*tau));
   %dsource_n=diag((-p.^2.-ui^2)./((p+n).^2.*tau));
   %DF=sparse([Mk(2:end-1, 2:end-1)*dsource_p, Mk(2:end-1, 2:end-1)*dsource_n; Mk(2:end-1, 2:end-1)*dsource_p, Mk(2:end-1, 2:end-1)*dsource_n]);

   density=y(1:N*N_species);
   DF=M*compute_change_rates_jacobian3(density, reactions, index,x); %jacobian matrix of the chemistry(N-2)* N_species x (N-2)* N_species


   jac_chem=sparse((N_species+1)*N,(N_species+1)*N);

   jac_chem(1:N_species*N, 1:N_species*N)=DF;

   jacrate=-L+jac_chem;

endfunction
