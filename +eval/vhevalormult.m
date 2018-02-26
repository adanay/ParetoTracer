% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [vHx, vHcount, vHundef,... 
          Hx, Hcount, Hundef] = vhevalormult(vH, H, x, Hx, v, iswarning, forcehess, opts)
% Vectorized Hessian multiply function evaluation.
% vH = [v' * H1; v' * H2; ...; v' * Hnobj].
% If vH is empty, it will perform the multiplication.
% Each row of x is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true.
% forcehess will force one evaluation of the Hessian to be returned in Hx 
% even if vH was passed and evaluated.
% opts are the optimization options.
%
% The result is always of size (m x nobj x n) even if m = 1.
% Hx (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).

forcehess = forcehess && nargout > 3;

vHundef = false;

Hcount = 0;
Hundef = false;

if isempty(vH)
  if isempty(Hx)
    [Hx, Hcount, Hundef] = eval.heval(H, x, iswarning, opts);
  end
  [vHx, vHcount] = eval.vhtimes(Hx, v);
else
  [vHx, vHcount, vHundef] = eval.vheval(vH, x, v, iswarning, opts);
end

if forcehess && isempty(Hx)
  [Hx, Hcount, Hundef] = eval.heval(H, x, iswarning, opts);
end
end