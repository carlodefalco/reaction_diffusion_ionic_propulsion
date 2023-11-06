clc
clear
close all
addpath (canonicalize_file_name ("../../../data"));
addpath (canonicalize_file_name ("./funzioni chimica manuali"));
addpath (canonicalize_file_name ("../IPROP_ARGON/READ_REACTIONS"));
[r, idx] = read_reactions (file_in_loadpath ("balcon_et_al_argon_ionization.json"));
N_species=numfields(idx);
pretty_print_reactions (r);

%===============================================================================
% LOAD BIM PACKAGE FOR FE MATRIX ASSEMBLY
%===============================================================================
pkg load fpl bim msh
%% MESH GENERATION (1D)
L = 2e-3; %[m]
N = 200;
x = linspace (0, L, N)';

state0 = zeros (numfields (idx), 1);
state0(idx.("Ar"))   = 2.5e+25;%[m^-3]`
state0(idx.("e"))    = 1.01e+18;%[m^-3]
state0(idx.("Ar+"))  = 1.0e18;%[m^-3]
state0(idx.("Ar2+")) = 1.0e16;%[m^-3]
state0(idx.("Ar*"))  = 1e14;%[m^-3]


y0 = zeros(N*(N_species),1);
distribution=@(x) (1- ((x-.9*L)/(.1*L)).^ 2) .* (x>.8*L & x<L)
y0(N*(idx.("e")-1)  +1:N*idx.("e"))   = (state0(idx.("e"))/0.066543).*distribution(x);
y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))  = state0(idx.("Ar"));
y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = (state0(idx.("Ar2+"))/0.066543).*distribution(x);
y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+")) = (state0(idx.("Ar+"))/0.066543).*distribution(x);
y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*")) = (state0(idx.("Ar*"))/0.066543).*distribution(x);

##y0(N*(idx.("e")-1)  +1:N*idx.("e"))     = state0(idx.("e"));
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))    = state0(idx.("Ar"));
##y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = state0(idx.("Ar2+"));
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+"))   = state0(idx.("Ar+"));
##y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*"))   = state0(idx.("Ar*"));
##y0(N*(idx.("e")-1)  +1:N*idx.("e"))     = state0(idx.("e")).* (x>.8*L & x<L);
##y0(N*(idx.("Ar")-1) +1:N*idx.("Ar"))    = state0(idx.("Ar")).* (x>.8*L & x<L);
##y0(N*(idx.("Ar2+")-1)+1:N*idx.("Ar2+")) = state0(idx.("Ar2+")).* (x>.8*L & x<L);
##y0(N*(idx.("Ar+")-1)+1:N*idx.("Ar+"))   = state0(idx.("Ar+")).* (x>.8*L & x<L);
##y0(N*(idx.("Ar*")-1)+1:N*idx.("Ar*"))   = state0(idx.("Ar*")).* (x>.8*L & x<L);

keyboard
T0   = 0;
Tend = 1;

Mk=bim1a_reaction(x, 1, 1); %Mass matrix for each species (NxN)
M=kron(eye(N_species), Mk);  %Mass matrix for all species ((N-2)*N_species x (N-2)*N_species)

dstate0 = reshape(compute_change_rates3 (y0, r, idx,x), N*N_species,1);

keyboard


fun = @(t, y, ydot) M*ydot - M*reshape(compute_change_rates3 (y, r, idx,x), N_species*N,1);
function [jy, jydot] = jacfun (t, y, ydot,x, M, r, idx, N, N_species)
  persistent t0 = 0
  if t > t0
    disp(t)
    t0 = t;
  endif
  jy = - M*compute_change_rates_jacobian3 (y, r, idx,x);
  jydot = M*eye(N*N_species);
endfunction

o = odeset ('AbsTol', 1e-4, 'RelTol', 1e-4, 'Jacobian', @(t, y, ydot) jacfun (t, y, ydot,x, M, r, idx,N, N_species), 'MaxOrder', 1);


keyboard

[t, y] = ode15i (fun, [0 logspace(-10, log10 (Tend), 1000)], y0, dstate0,o);
n1 = y(:,1:N);
n2 = y(:,N+1:2*N);
n3 = y(:,2*N+1:3*N);
n4 = y(:,3*N+1:4*N) ;
n5 = y(:,4*N+1:5*N) ;
n6 = y(:,5*N+1:6*N) ;

figure()
for ii = 1 :10: numel(t)

 semilogy(x, n1(ii, :),'LineWidth', 1,'r')
 hold on
 %semilogy(x, n2(ii, :)*1e-12,'LineWidth', 1,'b')
 semilogy(x, n3(ii, :),'LineWidth', 1,'b')
 semilogy(x, n4(ii, :),'LineWidth', 1,'m')
 semilogy(x, n5(ii, :),'LineWidth', 1,'c')
  semilogy(x, n6(ii, :),'LineWidth', 1,'y')
 axis ([min(x) max(x) 1e2 2.5e20]);
 str=sprintf('time = %d s', t(ii));
 legend('e', 'Ar+', 'Ar*', 'Ar2+', 'h_nu')
 title(str)
 hold off
 drawnow
endfor


##figure
##for [val, key] = idx
##  warning ("off", "Octave:negative-data-log-axis")
##  loglog (t(2:end), x(:, val)(2:end), 'Linewidth', 2)
##  set (gca, 'fontsize', 24,
##       'ytick', logspace(-12, 20, 9),
##       'xtick', logspace(-15, 5, 6))
##  ylim ([1e-12 1e26])
##  hold all
##endfor
##
##legend (fieldnames (idx){:}, 'location', 'eastoutside')
##hold off
