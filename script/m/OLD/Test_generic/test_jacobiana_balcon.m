addpath (canonicalize_file_name ("../../data"));
[r, idx] = read_reactions (file_in_loadpath ("balcon_et_al_argon_ionization.json"));

pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1);
x0(idx.("Ar"))   = 2.5e+19;
x0(idx.("e"))    = 1.0e+6;
x0(idx.("Ar+"))  = 1.0e+6;
x0(idx.("Ar2+")) = 1.0e+3;
x0(idx.("Ar*"))  = 1.0e+10;

T0   = 0;
Tend = 1.0e-7;

  N = numel (x0);
  d_Rfi_d_sj=zeros(1,N);
  d_Rbi_d_sj=zeros(1,N);
  dstate_jac = zeros (N, N);
 the_reaction = r(4)';
for [value, key] = the_reaction.reactants
  d_Rfi_d_sj(idx.(key)) = the_reaction.rate_coeffs(1) * x0(idx.(key)) ^ (value-1)* value;
  for [value2, key2] = the_reaction.reactants
    if key2!=key;
      d_Rfi_d_sj(idx.(key))*=x0(idx.(key2)) ^ (value2);
    else
      d_Rfi_d_sj(idx.(key))= d_Rfi_d_sj(idx.(key));
    endif
  endfor

endfor
for [value, key] = the_reaction.products
  d_Rbi_d_sj(idx.(key)) = the_reaction.rate_coeffs(2) * x0(idx.(key)) ^ (value-1)* value;
  for [value2, key2] = the_reaction.products
    if idx.(key2)!=idx.(key);
       d_Rbi_d_sj(idx.(key))*=x0(idx.(key2))^(value2);
    else
       d_Rbi_d_sj(idx.(key))= d_Rbi_d_sj(idx.(key));
    endif
  endfor
endfor

d_Ri_d_sj = d_Rfi_d_sj- d_Rbi_d_sj;
for [value, key] = the_reaction.reactants
    dstate_jac (idx.(key),:) -= d_Ri_d_sj * value;
endfor
for [value, key] = the_reaction.products
    dstate_jac (idx.(key),:) += d_Ri_d_sj * value;
endfor

