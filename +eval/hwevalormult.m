% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hwx, Hwcount, Hwundef,...
          Hx, Hcount, Hundef] = hwevalormult(Hw, H, x, Hx, w, iswarning, forcehess, opts)
% Vectorized Hessian multiply function evaluation.
% Hw = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj.
% If Hw is empty, it will perform the multiplication.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true.
% forcehess will force one evaluation of the Hessian to be returned in Hx 
% even if Hw was passed and evaluated.
% opts are the optimization options.
%
% The result is always of size (m x n x n) even if m = 1.
% Hx (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).

forcehess = forcehess && nargout > 3;

Hwundef = false;
Hcount = 0;
Hundef = false;

if isempty(Hw)
  if isempty(Hx)
    [Hx, Hcount, Hundef] = eval.heval(H, x, iswarning, opts);
  end
  [Hwx, Hwcount] = eval.hwtimes(Hx, w);
else
  [Hwx, Hwcount, Hwundef] = eval.hweval(Hw, x, w, iswarning, opts);
end

if forcehess && isempty(Hx)
  [Hx, Hcount, Hundef] = eval.heval(H, x, iswarning, opts);
end
end