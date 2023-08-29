addpath (canonicalize_file_name ("../../data"));
[r, idx] = read_reactions (file_in_loadpath ("robertson_autocatalysis.json"));

pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1);
x0(idx.("A"))   = 2.5e+19;
x0(idx.("B"))    = 1.0e+9;
x0(idx.("C"))  = 1.0e+6;


T0   = 0;
Tend = 1.0e-7;

  N = numel (x0);
  d_Rfi_d_sj=zeros(1,N);
  dstate_jac = zeros (N, N);
 the_reaction = r(1)';
for [value, key] = the_reaction.reactants
  d_Rfi_d_sj(idx.(key)) = the_reaction.rate_coeffs(1) * x0(idx.(key)) ^ (value-1)* value;
  for [value2, key2] = the_reaction.reactants
    if idx.(key2)!=idx.(key);
      d_Rfi_d_sj(idx.(key))*=x0(idx.(key2)) ^ (value2);
    else
       d_Rfi_d_sj(idx.(key))= d_Rfi_d_sj(idx.(key));
    endif
  endfor

endfor

for [value, key] = the_reaction.reactants
    dstate_jac (idx.(key),:) -= d_Rfi_d_sj * value;
endfor
for [value, key] = the_reaction.products
    dstate_jac (idx.(key),:) += d_Rfi_d_sj * value
endfor

%analytical derivatives of the rate coefficients for the robertson problem
% d_Rf1_d_sA=0.04; d_Rf1_d_sB=0;    d_Rf1_d_sC=0;
%d_Rf2_d_sA=0;    d_Rf2_d_sB=6e+7s_B;  d_Rf2_d_sC=0;
%d_Rf3_d_sA=0;  d_Rf3_d_sB=1e+4s_C;   d_Rf3_d_sC=1e+4s_B;

