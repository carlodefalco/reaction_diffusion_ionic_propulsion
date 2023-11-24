function [rate]=assembly_rate_transport(y, x, reactions, index, valence, mobility,bc,N, q, epsilon, Vth, P, Mk,k_b, me,T_e,nodal)
  N_species=numfields(index);
  N_species=N_species-1;
  %===============================================================================
  % RATE
  %===============================================================================

  A=sparse(N_species*N, N_species*N ); %global transport matrix (N_species*N x N_species*N )
  phi=y(N_species*N+1 : end);

  for k=1:N_species
    A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1,0);% -valence(k)*phi/Vth(k));
  endfor
   %%BOUNDARY CONDITIONS
##  A_bc=A;
  %A_bc(index.("e")*N,index.("e")*N )=A_bc(index.("e")*N,index.("e")*N)+ 0.1*(k_b*T_e/(2*pi*me))^0.5;
##  A_bc((index.("Ar+")-1)*N+1 ,:)=0;
##  A_bc((index.("Ar+")-1)*N+1,(index.("Ar+")-1)*N+1)=1;
##  A_bc((index.("Ar2+")-1)*N+1,:)=0;
##  A_bc((index.("Ar2+")-1)*N+1,(index.("Ar2+")-1)*N+1)=1;
####  A_bc([1 N], [1 N])=[];
##  %%test only diffusion
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
##  P([1 end],:)=0;
##  P(1,1)=1;
##  P(end, end)=1;
  K=kron(valence, -q.*Mk)';
  %K([1,end],[])=0;

  L=sparse([A, zeros(N_species*N, N);
            K , P]);%left-hand side matrix
##  L_bc=L;
##  L_bc(nodal, :)=[];
##  L_bc(:, nodal)=[];
##
  f1=sparse(N*(N_species+1),1); %right hand side matrix
##
##  for ii=1:numel(bc)/2
##  col=L(:,1+N*(ii-1)).*bc(2*ii-1);
##  f1-=col;
##  col=L(:,N*ii).*bc(ii*2);
##  f1-=col;
##  endfor
##f1(nodal)=[];
##  f1=-A(:,1)*bc(index.("e"),1)-A(:,N)*bc(index.("Ar"),2);
##  f1([(index.("e")-1)*N+1,index.("e")*N])=bc(index.("e"),:);
##  f1([(index.("Ar")-1)*N+1,index.("Ar")*N])=bc(index.("Ar"),:);
##  f1([(index.("Ar+")-1)*N+1,index.("Ar+")*N])=bc(index.("Ar+"),:);
##  f1([(index.("Ar2+")-1)*N+1,index.("Ar2+")*N])=bc(index.("Ar2+"),:);
##  f1([(index.("Ar*")-1)*N+1,index.("Ar*")*N])=bc(index.("Ar*"),:);
##  b_phi=zeros(N,1);
##  b_phi(1)=V0t(1);
##  b_phi(end)=V0t(2);

  %f1(N_species*N+1:end)= b_phi;
  rate =-L*y+f1;

endfunction
