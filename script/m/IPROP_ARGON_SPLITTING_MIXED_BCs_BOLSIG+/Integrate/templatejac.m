function [dfdy, dfdydot] = templatejac (dy, dydot, t, y, ydot)
    dfdy = dy (t, y, ydot) ;
    if (nargout>1)
      dfdydot = dydot (t, y, ydot) ;
    endif
endfunction
