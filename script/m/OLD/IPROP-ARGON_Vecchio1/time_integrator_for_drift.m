function [t,y]=time_integrator_for_drift( T_vec, x, y0, reactions, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth)


 N_species=numfields(index);
%===============================================================================
% CONSTANT MATRICES
%===============================================================================
 Mk=bim1a_reaction (x, 1, 1); %Mass matrix for each species (NxN)
 P= bim1a_laplacian (x, epsilon, 1); %stiffness matrix for phi(NxN)
 M=zeros(N_species*(N-2), N_species*(N-2)); %Mass matrix for all species (N*N_speciesxN*N_species)
 for k=1:N_species
    M(1+(N-2)*(k-1):(N-2)*k, 1+(N-2)*(k-1):(N-2)*k)=Mk(2:end-1, 2:end-1);
 endfor
%===============================================================================
% REMOVE BOUNDARY VALUES FROM THE INITIAL STATE
%===============================================================================
 y0_bc=[];
 for k=1:N_species
  y0_bc = [y0_bc;y0(2+(k-1)*N: k*N-1)];
 endfor
%===============================================================================
% ARRANGE Y IN A MATRIX FORM TO DISTINGUISH EACH SPECIES
%===============================================================================
% n_old is a matrix ((N-2)xN_species)
%rearrange the yvector according to the different species
 n_old=[];
 for k=1:N_species
   n_old= [n_old, y0_bc(1+(N-2)*(k-1):(N-2)*k)];
 endfor
%===============================================================================
% INITIALIZE ()_OLD VARIABLES
%===============================================================================
 y_old=y0_bc;

%===============================================================================
% INITIALIZE OUTPUT VARIABLES
%===============================================================================
 t=zeros(size(T_vec)); %instants at which the solution will be computed
 y=zeros(numel(T_vec), numel(y0)); %this matrix is (N_time x N_species*N)
 t(1)=T_vec(1); %initial instant
 y(1,:)=y0; %initial conditions
 V0t=bc_phi(T_vec(1)) %boundary condition for phi at initial time


  for ii=1:numel(T_vec)-1
%===============================================================================
% INTEGRATION IN TIME
%===============================================================================
   dt=T_vec(ii+1)-T_vec(ii); %compute the time step
%===============================================================================
% COMPUTE PHI(n_old)
%===============================================================================

   phi_old = zeros (N, 1);
   phi_old([1 end])= V0t;
   b=zeros(N-2,1);
   for k=1:N_species
     var=valence(k)*n_old(:,k);
     b+=var;
   endfor
   f_phi_no_bc=q*Mk(2:end-1, 2:end-1)*b;
   f_phi = f_phi_no_bc - P(2:end-1,1).*V0t(1)  -  P(2:end-1,end).*V0t(2);
   phi_old(2:end-1) = P(2:end-1, 2:end-1) \ f_phi;
%===============================================================================
% COMPUTE TRANSPORT MATRICES FOR EACH SPECIES USING PHI(n_old)
%===============================================================================

   A=zeros(N_species*N, N_species*N ); %global transport matrix (N_species*N x N_species*N )
   for k=1:N_species
     A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, -valence(k)*phi_old/Vth(k));
   endfor

   A_bc=zeros(N_species*(N-2),N_species*(N-2)); %remove the boundary raws and colums N_species*(N-2) x N_species*(N-2)
   for k=1:N_species
      A_bc(1+(N-2)*(k-1):(N-2)*k,1+(N-2)*(k-1):(N-2)*k) = A(2+N*(k-1):N*k-1,2+N*(k-1):N*k-1);
   endfor

%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================
   chemistry= compute_change_rates2 (n_old, reactions, index, x(2:end-1)); %chemical source term ((N-2)x N_species)
   F_old= zeros((N-2)*N_species,1); %compute source term wiht boundary conditions
   for k=1:N_species
     F_old(1+(N-2)*(k-1):k*(N-2))  = Mk(2:end-1, 2:end-1)* chemistry (:,k) - A(2+N*(k-1):N*k-1, N*(k-1)+1).*bc(k,1) - A(2+N*(k-1):N*k-1,N*k).*bc(k,2);
   endfor

   DF=M*compute_change_rates_jacobian2( n_old, reactions, index,x(2:end-1)); %jacobian matrix of the chemistry(N-2)* N_species x (N-2)* N_species


   L=(M./dt+A_bc-DF); %left-hand side matrix
   f=(M./dt-DF)*y_old+ F_old; %right hand side matrix
   y_new=L\f;
%===============================================================================
% ASSIGN THE COMPUTED VALUES TO THE OTUPUT VARIABLES
%===============================================================================
   t(ii+1)=T_vec(ii+1);
   %add the boundary values
   y_with_bc=zeros(N_species*N,1);
   for k=1:N_species
     y_with_bc(1+(k-1)*N:k*N)=[bc(k,1); y_new(1+(k-1)*(N-2):k*(N-2)) ;bc(k,2)];
   endfor


   y(ii+1, :)=y_with_bc;
   %%upload the ()_old variables with the ()_new ones
   y_old=y_new;
   n_old=[];
   for k=1:N_species
    n_old= [n_old, y_old(1+(N-2)*(k-1):(N-2)*k)];
   endfor
   V0t=bc_phi(T_vec(ii+1));
    t(ii+1)

 endfor

endfunction
