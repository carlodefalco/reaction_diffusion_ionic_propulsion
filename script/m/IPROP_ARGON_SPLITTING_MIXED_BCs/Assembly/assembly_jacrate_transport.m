%function [jacrate]=assembly_jacrate_transport(t, y, x, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk,ui, tau)
function [jacrate]=assembly_jacrate_transport(t,y, x,reactions, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk)
  N_species=numfields(index);
  M=kron(eye(N_species), Mk);
  V0t=bc_phi(t);
  %y_bc=y;

 A=sparse(N_species*N, N_species*N ); %global transport matrix (N_species*N x N_species*N )
 %A_bc=sparse(N_species*(N-2),N_species*(N-2));%remove the boundary raws and colums N_species*(N-2) x N_species*(N-2)
 %phi=[V0t(1); y_bc((N-2)*(N_species)+1:end); V0t(2) ];
 phi=y(N_species*N+1:end );
 for k=1:N_species
   A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, -valence(k)*phi/Vth(k));
   %A_bc(1+(N-2)*(k-1):(N-2)*k,1+(N-2)*(k-1):(N-2)*k) = A(2+N*(k-1):N*k-1,2+N*(k-1):N*k-1);
 endfor
  A_bc=A;
  A_bc((index.("Ar+")-1)*N+1 ,:)=0;
  A_bc((index.("Ar+")-1)*N+1,(index.("Ar+")-1)*N+1)=1;
  A_bc((index.("Ar2+")-1)*N+1,:)=0;
  A_bc((index.("Ar2+")-1)*N+1,(index.("Ar2+")-1)*N+1)=1;
  A_bc(index.("e")*N,index.("e")*N )=A_bc(index.("e")*N,index.("e")*N)+ Mk(end, end)*0.6;
  P([1 end],:)=0;
  P(1,1)=1;
  P(end, end)=1;
%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================


  L=sparse([(A_bc), zeros(N_species*N, N);
              kron(valence, -q.*Mk)', P]);%left-hand side matrix

  jacrate=-L;

endfunction
