addpath (canonicalize_file_name ("../../data"));
[r, idx] = read_reactions (file_in_loadpath ("balcon_et_al_argon_ionization.json"));

pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1);
x0(idx.("Ar"))   = 2.5e+19;
x0(idx.("e"))    = 1.0e+6;
x0(idx.("Ar+"))  = 1.0e+6;
x0(idx.("Ar2+")) = 1.0e+3;
x0(idx.("Ar*"))  = 1.0e+10;

T0   = 0.;
Tend = 1.0e-7;

[t, x] = ode23s (@(t, x) compute_change_rates (x, r, idx), [T0 Tend], x0);

figure
for [val, key] = idx
  semilogy (t(2:end), x(:, val)(2:end))
  hold all
endfor
legend (fieldnames (idx){:}, 'location', 'eastoutside')
hold off
