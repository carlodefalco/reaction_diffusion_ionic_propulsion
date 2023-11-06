close all;
clear all;
pkg load fpl bim msh

function u1 = lbwe_step (u0, t, dt, rate, jacrate, mass)
  b  = mass * u0 + dt * (rate (u0, t) - jacrate (u0, t) * u0) ;
  A  = mass - dt * jacrate (u0, t);
  u1 = full(A) \ b;
endfunction

function u1 = trap_step (u0, t, dt, rate, jacrate, mass)
  b  = mass * u0 + dt * (.001 * rate (u0, t) + .999 * (rate (u0, t) - jacrate (u0, t) * u0)) ;
  A  = mass - .999 * dt * jacrate (u0, t);
  u1 = full(A) \ b;
endfunction

%%mass (u1 - u0) / dt = 1/2 (rate (u0, t) )

function u = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass)

  u(:, 1) = u0;
  dt = (tspan(2) - tspan(1)) / 100;

  %%h = waitbar (0);
  for ispan = 1 : numel (tspan) - 1
    t  = tspan (ispan);
    ii = 0;
    while (t < tspan(ispan + 1))
      if (t + dt > tspan(ispan + 1) - 100*eps)
	tnew = tspan(ispan + 1);
        dt = tnew - t;
	t  = tnew;
      endif
      t = t + dt

      unew   = trap_step (u0, t, dt, rate, jacrate, mass);

      utest  = trap_step (u0, t-dt/2, dt/2, rate, jacrate, mass);
      utest  = trap_step (utest, t, dt/2, rate, jacrate, mass);
      err = norm ((unew - utest) ./ mag, inf)

      if (err < tol(t))
        u0 = unew;
        dt = .6 * (tol(t) / err)^(1/2) * dt
      else
        t = t-dt;
	printf ('reject dt_old = %g, dt_new = %g\n', dt, dt/2)
        dt = dt/2
      endif

    endwhile
    u(:, ispan+1) = unew;
    N = 300;
    L = 1.5e-6;
    x = linspace (0, L, N).';
    plot (x, diff (u(N+1:2*N, [ispan ispan+1]), 1, 2),
	  x, diff (u(2*N+1: 3*N, [ispan ispan+1]), 1, 2),
	  x, mag(2*N+1: 3*N) * tol(t) .* [-1, 1])
    drawnow
  endfor
  %%close (h)
endfunction


## u0 = [1; 1];
## rate = @(u, t) [-3*u(1); -2*u(2)];
## jacrate = @(u, t) [-3, 0; 0, -2];
## mass = eye(2);
## tol  = 1e-5;
## mag = 1;
## tspan = linspace (0, 10, 100);
## u = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass)


N = 300;
L = 1.5e-6;
x = linspace (0, L, N).';
D = 1e22 .* (x < .5*L);
A = 1e22 .* (x >= .5*L);
ni = 1e16;
Vth = 26e-3;
n = (D(1)+(-D(1)/2+sqrt((D(1)/2)^2+ni^2))*(exp(-.1/Vth))) .* (x <= .45*L)+(-A(end)/2+sqrt((A(end)/2)^2+ni^2)) .* (x >  .45*L);
p = (A(end)+(-A(end)/2+sqrt((A(end)/2)^2+ni^2))) .* (x >  .55*L)+(-D(1)/2+sqrt((D(1)/2)^2+ni^2))*(exp(-.1/Vth)) .* (x <= .55*L);

P = bim1a_laplacian (x, 4*8.85e-12, 1);
M = bim1a_reaction (x, 1, 1);

q = 1.6e-19;


mun = 1e-1;
mup = 1e-1;
tau = 1e-1;

phi = zeros (size (x));
%%phi(2:end-1) = P(2:end-1, 2:end-1) \ (q*M*(p-n+D-A))(2:end-1);
n(end) = ni;
p(1) = ni;
keyboard
function r = ratefun (u, t, N, x, P, q, M, D, A, Vth, mun, mup, tau, ni)

  phi = u(1:N);
  n   = u(N+1:2*N);
  p   = u(2*N+1:3*N);

  r1 = phi;
  r1(2:end-1) = P(2:end-1, 2:end-1) * phi(2:end-1) - (q*M*(p-n+D-A))(2:end-1);

  Cn = bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);
  r2 = zeros(N,1);
  r2([1 end]) = n([1 end]) - [(-D(1)/2+sqrt((D(1)/2)^2+ni^2))*(exp(-.1/Vth))+D(1); (-A(end)/2+sqrt((A(end)/2)^2+ni^2))];
  r2(2:end-1) = Cn(2:end-1, :) * n -(M*((ni^2 - p.*n)./((p+n)*tau)))(2:end-1);

  Cp = bim1a_advection_diffusion (x, mup*Vth, 1, 1, -phi/Vth);
  r3 = zeros(N,1);
  r3([1 end]) = p([1 end]) - [(-D(1)/2+sqrt((D(1)/2)^2+ni^2))*(exp(-.1/Vth)); (-A(end)/2+sqrt((A(end)/2)^2+ni^2))+A(end)];
  r3(2:end-1) = Cp(2:end-1, :) * p - (M*((ni^2 - p.*n)./((p+n)*tau)))(2:end-1);

  r = [-r1; -r2; -r3];

endfunction

function J = jacratefun (u, t, N, x, P, q, M, D, A, Vth, mun, mup, tau, ni)

  phi = u(1:N);
  n   = u(N+1:2*N);
  p   = u(2*N+1:3*N);

  Cn = bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);
  Cp = bim1a_advection_diffusion (x, mup*Vth, 1, 1, -phi/Vth);

  J11 = sparse (N,N);
  J11(1,1) = J11(N,N) = 1;
  J11(2:end-1, :) =  P(2:end-1, :);
  J12 = sparse (N,N);
  J12(2:end-1, :) =  q*M(2:end-1, :);
  J13 = sparse (N,N);
  J13(2:end-1, :) = -q*M(2:end-1, :);

  J21 = sparse (N,N);
  J22 = sparse (N,N);
  J22(1,1) = J22(N,N) = 1;
  J22(2:end-1, :) = Cn(2:end-1, :) + (M .* diag((ni^2+p.^2)./((p+n).^2*tau)))(2:end-1, :);
  J23 = sparse (N,N);
  J23(2:end-1, :) =  (M .* diag(  (ni^2+p.^2)./((p+n).^2*tau)  )  )(2:end-1, :);

  J31 = sparse (N,N);
  J32 = sparse (N,N);
  J32(2:end-1, :) = (M .* diag(  (ni^2+p.^2)./((p+n).^2*tau)  ) )(2:end-1, :);
  J33 = sparse (N,N);
  J33(1,1) = J33(N,N) = 1;
  J33(2:end-1, :) =  Cp(2:end-1, :) + (M .* diag(  (ni^2+n.^2)./((p+n).^2*tau)  )  )(2:end-1, :);

  J = -[J11,J12,J13; J21,J22,J23; J31,J32,J33];

endfunction

u0 = [phi; n; p];
tspan = linspace (0, 2e-10, 10);
mag = 100*ones(3*N, 1);
mag(N+1:N*2) = 1e22;
mag(2*N+1:N*3) = 1e22;
tol = @(t) .001 * (1-((tspan (1) - t) / (tspan(1) - tspan(end)))^.4) + 1e-6   * ((tspan (1) - t) / (tspan(1) - tspan(end)))^.4;
rate = @(u, t) ratefun (u, t, N, x, P, q, M, D, A, Vth, mun, mup, tau, ni);
jacrate = @(u, t) jacratefun (u, t, N, x, P, q, M, D, A, Vth, mun, mup, tau, ni);
mass = kron ([0 0 0; 0 1 0; 0 0 1], M);


u = integrate_adaptive (u0, tspan, tol, mag, rate, jacrate, mass);
figure
plot (x, u(1:N, :))
title('DEFALCO CODE')
figure
semilogy (x, u(N+1:2*N, end))
hold on
semilogy (x, u(2*N+1:3*N, end))
legend('n', 'p')
title('DEFALCO CODE')


