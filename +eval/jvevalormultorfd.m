function [Jvx, Jvcount, Jvundef, Jvapprox,...
          Jx, Jcount, Jundef, Japprox,...
          fx, fcount] = jvevalormultorfd(Jv, J, f, x, fx, Jx, v, lb, ub, iswarning, forcejac, opts)
% Vectorized Jacobian multiply function evaluation.
% Jv = J * v.
% If Jv is empty, it will perform the multiplication.
% If J is also empty, it will use finite differences (FD).
% Each row of x is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and Jv(x) or 
% J(x) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Jacobian will be computed by FD (the objective function must be available).
% + iswarning = false implies that an error will be thrown.
% forcejac will force one evaluation of the Jacobian to be returned in Jx 
% even if Jv was passed and evaluated.
% opts are the optimization options.
%
% The result is always of size (m x nobj) even if m = 1.
% The result of Jv(x) is assumed to be (nobj x 1) if m = 1.
% fx (if entered) is assumed to be of size (m x obj) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcejac = forcejac && nargout > 3;

Jvcount = 0;
Jcount = 0;
fcount = 0;

Jvundef = false;
Jvapprox = false;

Jundef = false;
Japprox = false;

if isempty(Jv) && isempty(J) && isempty(Jx)
  jacv();
else
  [Jvx, Jvcount, Jvundef, Jx, Jcount, Jundef] = eval.jvevalormult(Jv, J, x, Jx, v, iswarning, forcejac, opts);
  if Jvundef || Jundef
    Jx = [];
    jacv();
  end
end

  function [] = jacv()
    Jvapprox = true;
    if forcejac
      if isempty(Jx)
        [Jx, ~, ~, Japprox, fx, fcount] = eval.jevalorfd([], f, x, fx, lb, ub, iswarning, opts);
      end
      [Jvx, Jvcount] = eval.jvevalormult([], [], x, Jx, v, iswarning, forcejac, opts);
    else
      try
        [Jvx, fx, fcount] = fd.jacobian(f, x, fx, lb, ub, v, opts);
        m = max(size(x, 1), size(v, 1)); % m is the number of individuals
        if m == 1
          Jvx = Jvx(:)'; % (1 x nobj)
        end
      catch
        if isempty(Jx)
          [Jx, ~, ~, Japprox, fx, fcount] = eval.jevalorfd([], f, x, fx, lb, ub, iswarning, opts);
        end
        [Jvx, Jvcount] = eval.jvevalormult([], [], x, Jx, v, iswarning, forcejac, opts);
      end
    end
  end
end