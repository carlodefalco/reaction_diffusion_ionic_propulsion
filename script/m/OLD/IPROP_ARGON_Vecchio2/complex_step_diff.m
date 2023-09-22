function [J, DJ] = complex_step_diff ( s,  f)
  h = 1e-100;
  n = numel (s);

  for ii = 1 : n
    ds = s;
    ds(ii) += 1i * h;
    J(:, ii) = imag (f (ds)) / h;
  endfor

  x=linspace(0, 2e-3, 81);
  Mk=bim1a_reaction(x, 1, 1 );
  DJ=zeros(474,474);
  for k=1:6
    DJ(1+(81-2)*(k-1): k*(81-2), 1+(81-2)*(k-1):k*(81-2)) = Mk(2:end-1, 2:end-1);
  endfor

endfunction

%!test
%! f = @(x) [x(1).^2+3*exp(-x(2)/7); log(x(1) + 1).*x(2)];
%! Jex = @(x) [2*x(1), -3/7*exp(-x(2)/7); 1./(x(1) + 1).*x(2), log(x(1) + 1)];
%! s = rand (2, 1);
%! Jnum =  complex_step_diff (s, f);
%! disp (Jex(s))
%! disp (Jnum)
%! disp (norm (Jnum-Jex(s)))
