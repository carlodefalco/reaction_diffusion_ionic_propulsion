function [t,y, res_abs, res_rel]=time_integrator_for_drift( T_vec, x, y0, reactions, index, valence, mobility,bc, bc_phi,N, q, epsilon, Vth, P, Mk,ui, tau)


 N_species=numfields(index);
%===============================================================================
% CONSTANT MATRICES
%===============================================================================
  M=kron(eye(N_species), Mk(2:end-1, 2:end-1));  %Mass matrix for all species ((N-2)*N_species x (N-2)*N_species)
%===============================================================================
% REMOVE BOUNDARY VALUES FROM THE INITIAL STATE
%===============================================================================
 y0_bc=y0;
 y0_bc([1:N:(N_species+1)*N, N:N:(N_species+1)*N])=[];
%===============================================================================
% INITIALIZE ()_OLD VARIABLES.
 y_old=y0_bc;

%===============================================================================
% INITIALIZE OUTPUT VARIABLES
%===============================================================================
 t=zeros(size(T_vec)); %instants at which the solution will be computed
 y=zeros(numel(T_vec), numel(y0)); %this matrix is (N_time x N_species*N)
 t(1)=T_vec(1); %initial instant
 y(1,:)=y0; %initial conditions
 V0t=bc_phi(T_vec(1)) %boundary condition for phi at initial time


 res_abs=[];
 res_rel=[];
 iter_time=0;
  for ii=1:numel(T_vec)-1
tic
%===============================================================================
% INTEGRATION IN TIME
%===============================================================================
   dt=T_vec(ii+1)-T_vec(ii) %compute the time step

%===============================================================================
% COMPUTE TRANSPORT MATRICES FOR EACH SPECIES USING PHI(n_old)
%===============================================================================

   A=sparse(N_species*N, N_species*N ); %global transport matrix (N_species*N x N_species*N )
   A_bc=sparse(N_species*(N-2),N_species*(N-2));%remove the boundary raws and colums N_species*(N-2) x N_species*(N-2)
   phi_old=[V0t(1); y_old((N-2)*(N_species)+1:end); V0t(2) ];
   for k=1:N_species
     A(1+N*(k-1): k*N,1+N*(k-1): k*N) = bim1a_advection_diffusion(x, mobility(k)*Vth(k), 1, 1, -valence(k)*phi_old/Vth(k));
     A_bc(1+(N-2)*(k-1):(N-2)*k,1+(N-2)*(k-1):(N-2)*k) = A(2+N*(k-1):N*k-1,2+N*(k-1):N*k-1);
   endfor

%===============================================================================
% COMPUTE THE SOURCE TERM
%===============================================================================

   chemistry= compute_change_rates2 (y_old(1:(N-2)*N_species), reactions, index, x(2:end-1)); %chemical source term ((N-2)x N_species)

   F_old= zeros((N-2)*N_species,1); %compute source term wiht boundary conditions
   for k=1:N_species
     F_old(1+(N-2)*(k-1):k*(N-2))  = Mk(2:end-1, 2:end-1)* chemistry (:,k) - A(2+N*(k-1):N*k-1, N*(k-1)+1).*bc(k,1) - A(2+N*(k-1):N*k-1,N*k).*bc(k,2);
   endfor


   DF=M*compute_change_rates_jacobian2( y_old(1:(N-2)*N_species), reactions, index,x(2:end-1)); %jacobian matrix of the chemistry(N-2)* N_species x (N-2)* N_species




   L=sparse([(M./dt+A_bc-DF), zeros(N_species*(N-2), N-2);
              -q*kron(valence', ones(N-2)), P(2:end-1, 2:end-1)]);%left-hand side matrix


   f=sparse([(M./dt-DF)*y_old(1:(N_species*(N-2)))+ F_old;-P(2:end-1,1).*V0t(1)-P(2:end-1,end).*V0t(2)]); %right hand side matrix -DF


   y_new=L\f;

%===============================================================================
% ASSIGN THE COMPUTED VALUES TO THE OTUPUT VARIABLES
%===============================================================================
   t(ii+1)=T_vec(ii+1);
   %add the boundary values
   y_with_bc=zeros((N_species+1)*N,1);
   for k=1:N_species
     y_with_bc(1+(k-1)*N:k*N)=[bc(k,1); y_new(1+(k-1)*(N-2):k*(N-2)) ;bc(k,2)];
   endfor

   y_with_bc(N_species*N+1:end)=[V0t(1); y_new((N-2)*N_species+1:end); V0t(2)];


   y(ii+1, :)=y_with_bc;
   %%upload the ()_old variables with the ()_new ones
   y_old=y_new;
   V0t=bc_phi(T_vec(ii+1));

   res_abs_old=norm(f-L*y_new);
   res_rel_old=res_abs_old/norm(f);
   iter_time+=1

   res_abs=[res_abs;  res_abs_old];
   res_rel=[res_rel; res_rel_old];

toc

 endfor

endfunction
