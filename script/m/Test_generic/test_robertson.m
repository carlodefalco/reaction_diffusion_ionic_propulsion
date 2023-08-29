##this script is a test for the computation of the reaction rates for the ROBERTSON problem in a inmplicit formulation
##Two functions are adopted to compute the system of equations and the Jacobian.
##The DAE system is solved with ode15i

addpath (canonicalize_file_name ("../../data"));
addpath (canonicalize_file_name ("../"));
[r, idx] = read_reactions (file_in_loadpath ("robertson_autocatalysis.json"));

pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1); %by defaults the photon initial density is zero.
x0(idx.("A"))   = 1;
x0(idx.("B"))    = 0;
x0(idx.("C"))  = 0;



x0dot=zeros(numfields (idx), 1);

T0   = 0;
Tend = 1.0e-7;

eqs = @(t, x, xdot) compute_change_rates (x, xdot, r, idx);

options = odeset('RelTol', 10.0^(-7), 'AbsTol', 10.0^(-7), 'Jacobian', @compute_change_rates_jacobian);
[t, y] = ode15i ( eqs, [0 4*logspace(-6,6)], x0, x0dot);

figure

y(:,2) = 1e4*y(:,2);
semilogx(t,y)
ylabel('1e4 * y(:,2)')
title('Robertson DAE problem with a Conservation Law, solved by ODE15I')
