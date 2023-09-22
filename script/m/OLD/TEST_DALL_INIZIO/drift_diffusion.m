## Copyright (C) 2023 Carlo de Falco
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% Solve :
%%   n'- div (mun Vth (grad n - grad phi/Vth n))=n_chem
%%   p'- div (mun Vth (grad p + grad phi/Vth p))=p_chem
%%  -div (epsilon grad (phi)) = q * (p - n)
%%   n(0) = n(L) = p(0) = p(L) = 0
%%   phi(0) = 0 phi(L) = V0(t)

clc
clear

pkg load fpl bim msh
addpath (canonicalize_file_name ("../../../data"));
addpath (canonicalize_file_name ("../"));
addpath (canonicalize_file_name ("../Test_onlychemistry"));
[r, idx] = read_reactions (file_in_loadpath ("simple_reaction.json"));
elements=fieldnames(idx);
pretty_print_reactions (r);

L = 2e-3;
N = 81;
x = linspace (0, L, N).';
T0 = 0.999e-4;
T1 = 1.02e-4;
T_vec=linspace(T0, T1, 1e2);

function y0dot = initial_cond (t, y, x, N, V0, q, epsilon, Vth, mun, mup, reactions, index)
  V0t = V0(t);

  n   = y(1:N-2);
  p   = y((N-2)+1:2*(N-2));
  c   = y(2*(N-2)+1:3*(N-2));

  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  b   = q*(M(2:end-1,2:end-1)*(p - n + 0*c));

  phi = A00(2:end-1, 2:end-1) \ (b-A00(2:end-1, end) * V0t(2));
  phi=[ V0t(1); phi; V0t(2)];

  An =  bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);
  dn = An(2:end-1, 2:end-1) * n;


  Ap = bim1a_advection_diffusion (x, mup*Vth, 1, 1, phi/Vth);
  dp = Ap(2:end-1, 2:end-1) * p;

  Ac = bim1a_advection_diffusion (x, 1, 1, 1e-3, 0);
  dc = Ac(2:end-1, 2:end-1) * c;

  species=[n p c];
  chemistry=[];
  for ii=1:N-2
    chemistry=[chemistry; compute_change_rates(species(ii,:),reactions,index)] ;
  endfor
  f_1_no_bc   = chemistry (:,1);
  f_2_no_bc   = chemistry (:,2);
  f_3_no_bc   = chemistry (:,3);

  f_1   = M(2:end-1, 2:end-1)*f_1_no_bc - An(2:end-1,1)*0 - An(2:end-1, end)*0;
  f_2   = M(2:end-1, 2:end-1)*f_2_no_bc - Ap(2:end-1,1)*0 - Ap(2:end-1, end)*0;
  f_3   = M(2:end-1, 2:end-1)*f_3_no_bc - Ac(2:end-1,1)*1e17 - Ac(2:end-1, end)*1e17;
  y0dot=[M(2:end-1, 2:end-1)\(-dn+f_1);
          M(2:end-1, 2:end-1)\(-dp+f_2);
          M(2:end-1, 2:end-1)\(-dc+f_3)
          ];

endfunction

function eqs = odefun (t, y, ydot, x, N, V0, q, epsilon, Vth, mun, mup, reactions, index)

  V0t = V0(t);

  n   = y(1:N-2);
  p   = y((N-2)+1:2*(N-2));
  c   = y(2*(N-2)+1:3*(N-2));

  ndot   = ydot(1:N-2);
  pdot   = ydot((N-2)+1:2*(N-2));
  cdot   = ydot(2*(N-2)+1:3*(N-2));

  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  b   = q*(M(2:end-1,2:end-1)*(p - n));

  phi = A00(2:end-1, 2:end-1) \ (b-A00(2:end-1, end) * V0t(2));
  phi=[ V0t(1); phi; V0t(2)];

  An =  bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);
  dn = An(2:end-1, 2:end-1) * n;


  Ap = bim1a_advection_diffusion (x, mup*Vth, 1, 1, phi/Vth);
  dp = Ap(2:end-1, 2:end-1) * p;

  Ac = bim1a_advection_diffusion (x, 1, 1, 1e-3, 0);
  dc = Ac(2:end-1, 2:end-1) * c;
  species=[n p c];
  chemistry=[];
  for ii=1:N-2
    chemistry=[chemistry; compute_change_rates(species(ii,:),reactions,index)] ;
  endfor
  f_1_no_bc   = chemistry (:,1);
  f_2_no_bc   = chemistry (:,2);
  f_3_no_bc   = chemistry (:,3);
  f_1   = M(2:end-1, 2:end-1)*f_1_no_bc - An(2:end-1,1)*0 - An(2:end-1, end)*0;
  f_2   = M(2:end-1, 2:end-1)*f_2_no_bc - Ap(2:end-1,1)*0 - Ap(2:end-1, end)*0;
  f_3   = M(2:end-1, 2:end-1)*f_3_no_bc - Ac(2:end-1,1)*1e17 - Ac(2:end-1, end)*1e17;

  eqs=[M(2:end-1, 2:end-1)*ndot + dn-f_1;
       M(2:end-1, 2:end-1)*pdot + dp-f_2;
       M(2:end-1, 2:end-1)*cdot + dc-f_3

  ];

endfunction
function eqs = fun_for_complex (t, y, x, N, V0, q, epsilon, Vth, mun, mup, reactions, index)

  V0t = V0(t);

  n   = y(1:N-2);
  p   = y((N-2)+1:2*(N-2));
  c   = y(2*(N-2)+1:3*(N-2));

  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  b   = q*(M(2:end-1,2:end-1)*(p - n));

  phi = A00(2:end-1, 2:end-1) \ (b-A00(2:end-1, end) * V0t(2));
  phi=[ V0t(1); phi; V0t(2)];

  An =  bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);
  dn = An(2:end-1, 2:end-1) * n;


  Ap = bim1a_advection_diffusion (x, mup*Vth, 1, 1, phi/Vth);
  dp = Ap(2:end-1, 2:end-1) * p;

  Ac = bim1a_advection_diffusion (x, 1, 1, 1e-3, 0);
  dc = Ac(2:end-1, 2:end-1) * c;
  species=[n p c];
  chemistry=[];
  for ii=1:N-2
    chemistry=[chemistry; compute_change_rates(species(ii,:),reactions,index)] ;
  endfor
  f_1_no_bc   = chemistry (:,1);
  f_2_no_bc   = chemistry (:,2);
  f_3_no_bc   = chemistry (:,3);
  f_1   = M(2:end-1, 2:end-1)*f_1_no_bc - An(2:end-1,1)*0 - An(2:end-1, end)*0;
  f_2   = M(2:end-1, 2:end-1)*f_2_no_bc - Ap(2:end-1,1)*0 - Ap(2:end-1, end)*0;
  f_3   = M(2:end-1, 2:end-1)*f_3_no_bc - Ac(2:end-1,1)*1e17 - Ac(2:end-1, end)*1e17;

  eqs=[ An(2:end-1, 2:end-1) * n-f_1;
        Ap(2:end-1, 2:end-1) * p-f_2;
        Ap(2:end-1, 2:end-1) * p-f_3
       ];

endfunction

function [J,DJ] = jacobian (t, y, ydot, x, N, V0, q, epsilon, Vth, mun, mup, reactions, index)

  V0t = V0(t);

  n   = y(1:N-2);
  p   = y((N-2)+1:2*(N-2));
  c   = y(2*(N-2)+1:3*(N-2));

  ndot   = ydot(1:N-2);
  pdot   = ydot((N-2)+1:2*(N-2));
  cdot   = y(2*(N-2)+1:3*(N-2));


  A00 = bim1a_laplacian (x, epsilon, 1);
  M   = bim1a_reaction (x, 1, 1);
  b   = q*(M(2:end-1,2:end-1)*(p - n));

  phi = A00(2:end-1, 2:end-1) \ (b-A00(2:end-1, end) * V0t(2));
  phi=[ V0t(1); phi; V0t(2)];

  J = sparse (2*(N-2), 2*(N-2));
  DJ = sparse (2*(N-2), 2*(N-2));

  An = bim1a_advection_diffusion (x, mun*Vth, 1, 1, phi/Vth);
  Ap = bim1a_advection_diffusion (x, mup*Vth, 1, 1, phi/Vth);
  Ac = bim1a_advection_diffusion (x, 1, 1, 1e-3, 0);
  species=[n p c];
  df1_dn1=[];
  df1_dn2=[];
  df1_dn3=[];

  df2_dn1=[];
  df2_dn2=[];
  df2_dn3=[];

  df3_dn1=[];
  df3_dn2=[];
  df3_dn3=[];
   for ii=1:N-2
    f = @( y) compute_change_rates(y, reactions, index);
    J_N = complex_step_diff (species(ii,:), f);
    %J_N=compute_change_rates_jacobian( species(ii,:),reactions, index);


    df1_dn1=[df1_dn1; J_N(1,1)];
    df1_dn2=[df1_dn2; J_N(1,2)];
    df1_dn3=[df1_dn3; J_N(1,3)];

    df2_dn1=[df2_dn1; J_N(2,1)];
    df2_dn2=[df2_dn2; J_N(2,2)];
    df2_dn3=[df2_dn3; J_N(2,3)];

    df3_dn1=[df3_dn1; J_N(3,1)];
    df3_dn2=[df3_dn2; J_N(3,2)];
    df3_dn3=[df3_dn3; J_N(3,3)];

    endfor

  J = [An(2:end-1, 2:end-1)- diag(M(2:end-1,2:end-1)*df1_dn1),- diag(M(2:end-1,2:end-1)*df1_dn2),- diag(M(2:end-1,2:end-1)*df1_dn3);
        - diag(M(2:end-1,2:end-1)*df2_dn1),Ap(2:end-1, 2:end-1)-diag(M(2:end-1,2:end-1)*df2_dn2),- diag(M(2:end-1,2:end-1)*df2_dn3);
         - diag(M(2:end-1,2:end-1)*df3_dn1),-diag(M(2:end-1,2:end-1)*df3_dn2),Ac(2:end-1, 2:end-1)- diag(M(2:end-1,2:end-1)*df3_dn3)
      ];


  DJ(1:N-2, 1:N-2) = M(2:end-1, 2:end-1);
  DJ((N-2)+1:2*(N-2), (N-2)+1:2*(N-2)) = M(2:end-1, 2:end-1);
  DJ(2*(N-2)+1:3*(N-2), 2*(N-2)+1:3*(N-2)) = M(2:end-1, 2:end-1);

endfunction
function [J2,DJ] = jacobian2 (t, y, ydot, x, N, V0, q, epsilon, Vth, mun, mup, reactions, index)
  M=bim1a_reaction(x,1,1);
  T=t;
  f2=@(ydot) odefun(T, y, ydot, x, N, V0, q, epsilon, Vth, mun, mup, reactions, index);
  J2=complex_step_diff(ydot, f2);

  DJ(1:N-2, 1:N-2) = M(2:end-1, 2:end-1);
  DJ((N-2)+1:2*(N-2), (N-2)+1:2*(N-2)) = M(2:end-1, 2:end-1);
  DJ(2*(N-2)+1:3*(N-2), 2*(N-2)+1:3*(N-2)) = M(2:end-1, 2:end-1);

endfunction


V0 = @(t) 5000*[zeros(size(t)); ((t-1e-4)*1e6 .* ((t>=1e-4)&(t<=1.01e-4)) + 1.0 .* (t>1.01e-4))];
q  = 1.6e-19;
epsilon = 8.8e-12;
mun = 1e-3;
mup = 1e-3;
Vth = 26e-3;

n0 = x .* (L - x) * 1e17 / (L/2)^2;
p0 = zeros(size(x));
c0= ones(N,1)*1e17;

y0 = [n0(2:end-1); p0(2:end-1);c0(2:end-1)];
y0dot= initial_cond (T0, y0, x, N, V0, q, epsilon, Vth, mun, mup,r,idx);
 [J2,DJ] = jacobian2 (T0, y0, y0dot, x, N, V0, q, epsilon, Vth, mun, mup, r, idx);
 [J,DJ] = jacobian (T0, y0, y0dot, x, N, V0, q, epsilon, Vth, mun, mup, r, idx);

o = odeset ('RelTol', 1e-7, 'AbsTol', 1e-8, 'InitialStep', 1e-25,...
'Jacobian', @(t, y, ydot) jacobian (t, y,ydot, x, N, V0, q, epsilon, Vth, mun, mup, r, idx));
fun=@(t, y, ydot) odefun (t, y, ydot, x, N, V0, q, epsilon, Vth, mun, mup, r, idx);
[t, y] = ode15i (fun,T_vec, y0, y0dot, o);

n = [zeros(size(t)),y(:, 1:N-2),zeros(size(t))];
p = [zeros(size(t)),y(:, (N-2)+(1:N-2)),zeros(size(t))];
c = [ones(size(t))*1e17,y(:, 2*(N-2)+(1:N-2)),ones(size(t))*1e17];
V = V0(t')(2, :);
for ii = 1 :1:numel (t)
  subplot(1, 2, 1)
  plot ( x, p(ii, :), x, c(ii, :))%x,n(ii, :)
  title (sprintf ("%g", t(ii)));
  legend ('n', 'p', 'c');
  axis ([min(x) max(x) 0 max(max(y0))]);
  subplot(1, 2, 2)
  plot (t, V, t(ii), V(ii), 'ro')
  title ('V')
  xlabel ('t [s]')
  ylabel ('V [V]')
  drawnow
endfor
