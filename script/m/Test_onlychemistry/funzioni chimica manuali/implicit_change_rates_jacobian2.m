function [dstate_jac, ddstate_jac] = implicit_change_rates_jacobian2(t, state, dstate, reactions, index)
  s_A=state(index.("A"));
  s_B=state(index.("B"));
  s_C=state(index.("C"));

  s_A_dot = dstate(index.("A"));
  s_B_dot = dstate(index.("A"));
  s_C_dot = dstate(index.("A"));

  dstate_jac=[ -reactions(1).rate_coeffs(1), reactions(3).rate_coeffs(1)*s_C,                                    reactions(3).rate_coeffs(1)*s_B;
               -reactions(1).rate_coeffs(1), -(reactions(3).rate_coeffs(1)*s_C+2*reactions(2).rate_coeffs(1)*s_B),-reactions(3).rate_coeffs(1)*s_B;
               0                           , -2*reactions(2).rate_coeffs(1)*s_B                                   , 0
              ]
    ddstate_jac= sparse(eye(N));
    endfunction
