% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Jx, Jcount, Jundef, Japprox, fx, fcount] = jevalorfd(J, f, x, fx, lb, ub, iswarning, opts)
% Vectorized Jacobian function evaluation.
% If J is empty, finite differences (FD) are used.
% Each row of x is considered an individual to be evaluated.
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and J(x) is 
% undef, then
% + iswarning = true implies that a warning will be displayed and the
% Jacobian will be computed by FD (the objective function must be available).
% + iswarning = false implies that an error will be thrown.
% opts are the FD and optimization options.
%
% The result is always of size (m x nobj x n) even if m = 1.
% The result of J(x) is assumed to be (nobj x n) if m = 1.
% fx (if entered) is assumed to be of size (m x obj) even if m = 1 
% (after being calculated using one of the vec eval functions).

Jcount = 0;
fcount = 0;
Jundef = false;
Japprox = false;

if isempty(J)
  jac();
else
  [Jx, Jcount, Jundef] = eval.jeval(J, x, iswarning, opts);
  if Jundef
    jac();
  end
end

  function [] = jac()
    Japprox = true;
    [Jx, fx, fcount] = fd.jacobian(f, x, fx, lb, ub, [], opts);
    m = size(x, 1); % m is the number of individuals
    if m == 1
      Jx = shiftdim(Jx, -1); % (1 x nobj x n)
    end
  end
end