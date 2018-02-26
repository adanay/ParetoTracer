% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [wJx, wJcount, wJundef,...
          Jx, Jcount, Jundef] = wjevalormult(wJ, J, x, Jx, w, iswarning, forcejac, opts)
% Vectorized Jacobian multiply function evaluation.
% wJ = w' * J.
% If wJ is empty, it will perform the multiplication.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true.
% forcejac will force one evaluation of the Jacobian to be returned in Jx 
% even if wJ was passed and evaluated.
% opts are the optimization options.
%
% The result is always of size (m x n) even if m = 1.
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).

forcejac = forcejac && nargout > 3;

Jcount = 0;
wJundef = 0;
Jundef = false;

if isempty(wJ)
  if isempty(Jx)
    [Jx, Jcount, Jundef] = eval.jeval(J, x, iswarning, opts);
  end
  [wJx, wJcount] = eval.wjtimes(Jx, w);
else
  [wJx, wJcount, wJundef] = eval.wjeval(wJ, x, w, iswarning, opts);
end

if forcejac && isempty(Jx)
  [Jx, Jcount, Jundef] = eval.jeval(J, x, iswarning, opts);
end
end