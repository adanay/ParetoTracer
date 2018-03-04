function [Jvx, Jvcount, Jvundef,...
          Jx, Jcount, Jundef] = jvevalormult(Jv, J, x, Jx, v, iswarning, forcejac, opts)
% Vectorized Jacobian multiply function evaluation.
% Jv = J * v.
% If Jv is empty, it will perform the multiplication.
% Each row of x is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true.
% forcejac will force one evaluation of the Jacobian to be returned in Jx 
% even if Jv was passed and evaluated.
% opts are the optimization options.
%
% The result is always of size (m x nobj) even if m = 1.
% The result of Jv(x) is assumed to be (nobj x 1) if m = 1.
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcejac = forcejac && nargout > 3;

Jcount = 0;
Jvundef = false;
Jundef = false;

if isempty(Jv)
  if isempty(Jx)
    [Jx, Jcount, Jundef] = eval.jeval(J, x, iswarning, opts);
  end
  [Jvx, Jvcount] = eval.jvtimes(Jx, v);
else
  [Jvx, Jvcount, Jvundef] = eval.jveval(Jv, x, v, iswarning, opts);
end

if forcejac && isempty(Jx)
  [Jx, Jcount, Jundef] = eval.jeval(J, x, iswarning, opts);
end
end