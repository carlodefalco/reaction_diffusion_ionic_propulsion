function [dstate_jac, ddstate_jac] = implicit_change_rates_jacobian_balcon(t,state, dstate, reactions, index)
  N_species=numfields(index);

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

  dR1_f_s_e  = reactions(1).rate_coeffs(1)*s_Ar;
  dR1_f_s_Ar = reactions(1).rate_coeffs(1)*s_e;

  dR2_f_s_e = reactions(2).rate_coeffs(1)*s_Ar;
  dR2_f_s_Ar = reactions(2).rate_coeffs(1)*s_e;

  dR3_f_s_e  = reactions(3).rate_coeffs(1)*s_Ar_excited;
  dR3_f_s_Ar_excited = reactions(3).rate_coeffs(1)*s_e;


  dR4_f_s_Ar_excited = 2 * reactions(4).rate_coeffs(1)*s_Ar_excited;

  dR5_f_s_Ar  = 2*reactions(5).rate_coeffs(1)*s_Ar*s_Ar_plus;
  dR5_f_s_Ar_plus = reactions(5).rate_coeffs(1)*s_Ar^2;

  dR6_f_s_e  = reactions(6).rate_coeffs(1)*s_Ar2_plus;
  dR6_f_s_Ar2_plus = reactions(6).rate_coeffs(1)*s_e;

  dR7_f_s_Ar_excited  = reactions(7).rate_coeffs(1);


  dR8_f_s_e  = reactions(8).rate_coeffs(1)*s_Ar;
  dR8_f_s_Ar = reactions(8).rate_coeffs(1)*s_e;


  sum_dRi_e_e = dR1_f_s_e + dR3_f_s_e -dR6_f_s_e;
  sum_dRi_e_Ar= dR1_f_s_Ar ;
  sum_dRi_e_Arplus=0;
  sum_dRi_e_Ar2plus=-dR6_f_s_Ar2_plus;
  sum_dRi_e_Arexcited= dR3_f_s_Ar_excited+dR4_f_s_Ar_excited;
  sum_dRi_e_hnu=0;

  sum_dRi_Ar_e  = -dR1_f_s_e - dR2_f_s_e + dR6_f_s_e ;
  sum_dRi_Ar_Ar = - dR1_f_s_Ar - dR2_f_s_Ar -dR5_f_s_Ar ;
  sum_dRi_Ar_Arplus =  -dR5_f_s_Ar_plus ;
  sum_dRi_Ar_Ar2plus =  + dR6_f_s_Ar2_plus ;
  sum_dRi_Ar_Arexcited =   dR4_f_s_Ar_excited+ dR7_f_s_Ar_excited;
  sum_dRi_Ar_hnu = 0;

  sum_dRi_Ar_plus_e = dR1_f_s_e +dR3_f_s_e  ;
  sum_dRi_Ar_plus_Ar =dR1_f_s_Ar  - dR5_f_s_Ar;
  sum_dRi_Ar_plus_Arplus =  - dR5_f_s_Ar_plus;
  sum_dRi_Ar_plus_Ar2plus=  0;
  sum_dRi_Ar_plus_Arexcited = + dR3_f_s_Ar_excited  +dR4_f_s_Ar_excited ;
  sum_dRi_Ar_plus_hnu = 0;

  sum_dRi_Ar2_plus_e =  - dR6_f_s_e ;
  sum_dRi_Ar2_plus_Ar = + dR5_f_s_Ar;
  sum_dRi_Ar2_plus_Arplus =  dR5_f_s_Ar_plus ;
  sum_dRi_Ar2_plus_Ar2plus = - dR6_f_s_Ar2_plus;
  sum_dRi_Ar2_plus_Arexcited =  0;
  sum_dRi_Ar2_plus_hnu =  0;

  sum_dRi_Ar_excited_e = + dR2_f_s_e - dR3_f_s_e  +dR6_f_s_e  ;
  sum_dRi_Ar_excited_Ar = +dR2_f_s_Ar   ;
  sum_dRi_Ar_excited_Arplus =  0;
  sum_dRi_Ar_excited_Ar2plus = + dR6_f_s_Ar2_plus ;
  sum_dRi_Ar_excited_Arexcited =  - dR3_f_s_Ar_excited- 2*dR4_f_s_Ar_excited  - dR7_f_s_Ar_excited ;
  sum_dRi_Ar_excited_hnu = +0;

  sum_dRi_hnu_e = 0;
  sum_dRi_hnu_Ar = 0;
  sum_dRi_hnu_Arplus = 0;
  sum_dRi_hnu_Ar2plus = 0;
  sum_dRi_hnu_Arexcited = dR7_f_s_Ar_excited ;
  sum_dRi_hnu_hnu = 0;

  sum_dRi_e           =[sum_dRi_e_e         ,sum_dRi_e_Ar         ,sum_dRi_e_Arplus         ,sum_dRi_e_Arexcited           ,sum_dRi_e_Ar2plus         ,         sum_dRi_e_hnu];
  sum_dRi_Ar          =[sum_dRi_Ar_e        ,sum_dRi_Ar_Ar        ,sum_dRi_Ar_Arplus        ,sum_dRi_Ar_Arexcited          ,sum_dRi_Ar_Ar2plus        ,        sum_dRi_Ar_hnu];
  sum_dRi_Ar_plus     =[sum_dRi_Ar_plus_e   ,sum_dRi_Ar_plus_Ar   ,sum_dRi_Ar_plus_Arplus   ,sum_dRi_Ar_plus_Arexcited     ,sum_dRi_Ar_plus_Ar2plus   ,   sum_dRi_Ar_plus_hnu];
  sum_dRi_Ar2_plus    =[sum_dRi_Ar2_plus_e  ,sum_dRi_Ar2_plus_Ar  ,sum_dRi_Ar2_plus_Arplus  ,sum_dRi_Ar2_plus_Arexcited    ,sum_dRi_Ar2_plus_Ar2plus  ,  sum_dRi_Ar2_plus_hnu];
  sum_dRi_Ar_excited  =[sum_dRi_Ar_excited_e,sum_dRi_Ar_excited_Ar,sum_dRi_Ar_excited_Arplus,sum_dRi_Ar_excited_Arexcited  ,sum_dRi_Ar_excited_Ar2plus,sum_dRi_Ar_excited_hnu];
  sum_dRi_hnu         =[sum_dRi_hnu_e       ,sum_dRi_hnu_Ar       ,sum_dRi_hnu_Arplus       ,sum_dRi_hnu_Arexcited         ,sum_dRi_hnu_Ar2plus       ,       sum_dRi_hnu_hnu];

  dstate_jac=-[sum_dRi_e;
               sum_dRi_Ar;
               sum_dRi_Ar_plus;
               sum_dRi_Ar_excited;
               sum_dRi_Ar2_plus;
               sum_dRi_hnu
              ];
  ddstate_jac= sparse(eye(N_species));
endfunction
