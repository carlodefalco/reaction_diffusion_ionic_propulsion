function [t,y]=time_integrator( T_vec, state0, reactions, index)
 %system assembly
 N_species=numfields(index);
 M=eye(N_species);
 dstate_old=compute_change_rates(state0, reactions, index);
 ddstate=compute_change_rates_jacobian(state0, reactions, index);
 state_old=state0;
 %%initialization of the output array
 t=zeros(size(T_vec));
 y=zeros(numel(T_vec), numel(state0));
 t(1)=T_vec(1);
 y(1,:)=state0;
%%time integration
 for ii=1:numel(T_vec)-1
   dt=T_vec(ii+1)-T_vec(ii);
   L=(M./dt-ddstate);
   b=(M./dt-ddstate)*state_old+dstate_old;
   state_new=L\b;
   %% assaign the computed value to the output
   t(ii+1)=T_vec(ii+1);
   y(ii+1, :)=state_new;
   %%upload the ()_old variables with the ()_new ones
   state_old=state_new;
   dstate_old=compute_change_rates(state_old, reactions, index);
   ddstate=compute_change_rates_jacobian(state_old, reactions, index);

 endfor
endfunction
