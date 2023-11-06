##this script is a test for the computation of the reaction rates for the ROBERTSON problem in a inmplicit formulation
##Two functions are adopted to compute the system of equations and the Jacobian.
##The DAE system is solved with ode15i


addpath (canonicalize_file_name ("../../../data"));
addpath (canonicalize_file_name ("./funzioni chimica manuali"));
addpath (canonicalize_file_name ("../"));
[r, idx] = read_reactions (file_in_loadpath ("robertson_autocatalysis.json"));
pretty_print_reactions (r);

x0 = zeros (numfields (idx), 1); %by defaults the photon initial density is zero.
x0(idx.("A"))   = 1;
x0(idx.("B"))    = 0.5;
x0(idx.("C"))  = 0.5;

x0dot_hat=[-0.04,0.04,0]';
T= [0 4*logspace(-6,6, 1e3)];
T_hat=T;%./t_bar;


x0dot=zeros(numfields (idx), 1);

##eqs =compute_change_rates_implicit2 (x0_hat, x0dot_hat, r, idx);
##eqs3 =compute_change_rates_implicit (x0_hat, x0dot_hat, r, idx);
##eqs2 = @(t, x, xdot) compute_change_rates_implicit2 (x, xdot, r, idx);
##
##options = odeset('RelTol', 10.0^(-7), 'AbsTol', 10.0^(-7), 'Jacobian', @implicit_change_rates_jacobian);
##[t, y] = ode15i ( eqs2, T_hat, x0, x0dot);
J=compute_change_rates_jacobian(x0, r, idx);
keyboard
[t3,y3]= time_integrator( T, x0, r, idx);

figure()
y2 = y;%.*x_bar;
y2(:,2) = y2(:,2)*1e4;
semilogx(t,y2)


hold on
y2 = y3;%.*x_bar;
y2(:,2) = y2(:,2)*1e4;
semilogx(t,y2)
ylabel('1e4 * y(:,2)')
legend('A-ode15i', 'B-ode15i','C-ode15i','A-manual', 'B-manual','C-manual')
title('Robertson DAE problem ')

