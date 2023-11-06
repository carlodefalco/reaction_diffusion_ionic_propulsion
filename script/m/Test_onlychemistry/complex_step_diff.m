function [J, DJ] = complex_step_diff ( s,  f)
  h = 1e-100;
  n = numel (s);

  for ii = 1 : n
    ds = s;
    ds(ii) += 1i * h;
    J(:, ii) = imag (f (ds)) / h;
  endfor
  DJ=eye(6);
endfunction

%!test
%! f = @(x) [x(1).^2+3*exp(-x(2)/7); log(x(1) + 1).*x(2)];
%! Jex = @(x) [2*x(1), -3/7*exp(-x(2)/7); 1./(x(1) + 1).*x(2), log(x(1) + 1)];
%! s = rand (2, 1);
%! Jnum =  complex_step_diff (s, f);
%! disp (Jex(s))
%! disp (Jnum)
%! disp (norm (Jnum-Jex(s)))
