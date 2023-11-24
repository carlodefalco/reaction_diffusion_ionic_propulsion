function [jacrate]=assembly_jacrate_transport(y, x,reactions, index, valence, mobility,bc,N, q, epsilon, Vth, P, Mk, k_b, me, T_e,nodal)
 N_species=numfields(index);
 N_species=N_species-1;
 M=kron(eye(N_species), Mk);
 %V0t=bc_phi;

 A=sparse(N_species*N, N_species*N ); %global transport matrix (N_species*N x N_species*N )
 phi=y(N_species*N+1:end );
 for k=1:N_species
   A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, 0);%-valence(k)*phi/Vth(k));
 endfor
 A_bc=A;
## A_bc((index.("Ar+")-1)*N+1 ,:)=0;
## A_bc((index.("Ar+")-1)*N+1,(index.("Ar+")-1)*N+1)=1;
## A_bc((index.("Ar2+")-1)*N+1,:)=0;
## A_bc((index.("Ar2+")-1)*N+1,(index.("Ar2+")-1)*N+1)=1;
## %A_bc(index.("e")*N,index.("e")*N )=A_bc(index.("e")*N,index.("e")*N)+  0.1*(k_b*T_e/(2*pi*me))^0.5;
##
##   %%test only diffusion
##  A_bc((index.("e")-1)*N+1 ,:)=0;
##  A_bc((index.("e")-1)*N+1,(index.("e")-1)*N+1)=1;
##  A_bc((index.("Ar")-1)*N+1 ,:)=0;
##  A_bc((index.("Ar")-1)*N+1,(index.("Ar")-1)*N+1)=1;
##  A_bc((index.("Ar*")-1)*N+1 ,:)=0;
##  A_bc((index.("Ar*")-1)*N+1,(index.("Ar*")-1)*N+1)=1;
##
##  A_bc((index.("e"))*N ,:)=0;
##  A_bc(index.("e")*N,index.("e")*N)=1;
##  A_bc((index.("Ar"))*N ,:)=0;
##  A_bc(index.("Ar")*N,index.("Ar")*N)=1;
##  A_bc((index.("Ar+"))*N ,:)=0;
##  A_bc(index.("Ar+")*N,index.("Ar+")*N)=1;
##  A_bc((index.("Ar2+"))*N ,:)=0;
##  A_bc(index.("Ar2+")*N,index.("Ar2+")*N)=1;
##  A_bc((index.("Ar*"))*N ,:)=0;
##  A_bc(index.("Ar*")*N,index.("Ar*")*N)=1;
## P([1 end],:)=0;
## P(1,1)=1;
## P(end, end)=1;
 K=kron(valence, -q.*Mk)';
## K([1 end],end)=0;
%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================
 L=sparse([(A_bc), zeros(N_species*N, N);
              K, P]);%left-hand side matrix
##  L_bc=L;
##  L_bc(nodal, :)=[];
##  L_bc(:, nodal)=[];
 jacrate=-L;

endfunction
