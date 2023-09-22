addpath (canonicalize_file_name ("../../data"));
[r, idx] = read_reactions (file_in_loadpath ("robertson_autocatalysis.json"));

pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1);
x0(idx.("A"))   = 1;
x0(idx.("B"))    = 0;
x0(idx.("C"))  = 0;


T0   = 0;
Tend = 1.0e-7;
J = @(t, x)  compute_change_rates_jacobian(x, r, idx);
J(0,x0)
%analytical derivatives of the rate coefficients for the robertson problem
% d_sA'_d_sA=-0.04; d_sA'_d_sB= 1.0e+10;    d_sA'_d_sC=1.0e+10;
%d_sB'_d_sA=0.04;    d_sB'_d_sB=-6,001e+13;  d_sB'_d_sC=-1.0e+10;
%d_sC'_d_sA=0;  d_sC'_d_sB=6.0e+13;   d_sC'_d_sC=0;

