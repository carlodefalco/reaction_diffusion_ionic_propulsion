function Res = compute_change_rates_implicit_balcon (state, dstate, reactions, index)

  s_e =state(index.("e"));
  s_Ar=state(index.("Ar"));
  s_Ar_plus =state(index.("Ar+"));
  s_Ar2_plus=state(index.("Ar2+"));
  s_Ar_excited=state(index.("Ar*"));
  s_hnu=state(index.("h_nu"));

  s_e_dot =dstate(index.("e"));
  s_Ar_dot=dstate(index.("Ar"));
  s_Ar_plus_dot =dstate(index.("Ar+"));
  s_Ar2_plus_dot=dstate(index.("Ar2+"));
  s_Ar_excited_dot=dstate(index.("Ar*"));
  s_hnu_dot=dstate(index.("h_nu"));

  R1_f = reactions(1).rate_coeffs(1)*s_e*s_Ar;
  R2_f = reactions(2).rate_coeffs(1)*s_e*s_Ar;
  R3_f = reactions(3).rate_coeffs(1)*s_e*s_Ar_excited;
  R4_f = reactions(4).rate_coeffs(1)*s_Ar_excited^2;
  R5_f = reactions(5).rate_coeffs(1)*s_Ar_plus*s_Ar^2;
  R6_f = reactions(6).rate_coeffs(1)*s_e*s_Ar2_plus;
  R7_f = reactions(7).rate_coeffs(1)*s_Ar_excited;
  R8_f = reactions(8).rate_coeffs(1)*s_e*s_Ar;

##  R1_b = 0;
##  R2_b = 0;
##  R3_b = 0;
##  R4_b = 0;
##  R5_b = 0;
##  R6_b = 0;
##  R7_b = 0;
##  R8_b = 0;


  sum_Ri_e = R1_f + R3_f + R4_f - R6_f;
  sum_Ri_Ar = - R1_f - R2_f + R4_f -R5_f + R6_f + R7_f;
  sum_Ri_Ar_plus = R1_f + R3_f + R4_f - R5_f;
  sum_Ri_Ar2_plus = + R5_f - R6_f;
  sum_Ri_Ar_excited = + R2_f - R3_f- 2*R4_f + R6_f - R7_f;
  sum_Ri_hnu = R7_f;

  Res=[ s_e_dot           - sum_Ri_e;
        s_Ar_dot          - sum_Ri_Ar;
        s_Ar_plus_dot     - sum_Ri_Ar_plus;
        s_Ar_excited_dot  - sum_Ri_Ar_excited;
        s_Ar2_plus_dot    - sum_Ri_Ar2_plus;
        s_hnu_dot         - sum_Ri_hnu ;
       ];





endfunction
