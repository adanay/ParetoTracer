% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [wJvx, wJvcount, wJvundef, wJvapprox,... 
          wJx, wJcount, wJundef, wJapprox,... 
          Jvx, Jvcount, Jvundef, Jvapprox,...
          Jx, Jcount, Jundef, Japprox,...
          wfx, wfcount, fx, fcount] =...
  wjvevalormultorfd(wJv, wJ, Jv, J, f, x, wfx, fx, wJx, Jvx, Jx, w, v, lb, ub, iswarning, forcejac, opts)
% Vectorized Jacobian multiply function evaluation.
% wJv = w' * J * v.
% If wJv is empty, it will perform the multiplication.
% If wJ, Jv, and J are also empty, it will use finite differences.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x, w, and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and wJv(x) or 
% wJ(x) or Jv(x) or J(x) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Jacobian will be computed by FD (the objective function must be available).
% + iswarning = false implies that an error will be thrown.
% forcejac will force one evaluation of the Jacobian to be returned in Jx 
% even if wJv, wJ, or Jv, were passed and evaluated.
% opts are the optimization options.
%
% The result is always of size (m x 1).
% fx (if entered) is assumed to be of size (m x obj) even if m = 1 
% (after being calculated using one of the vec eval functions).
% wfx (if entered) is assumed to be of size (m x 1) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% wJx (if entered) is assumed to be of size (m x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Jvx (if entered) is assumed to be of size (m x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).

forcejac = forcejac && nargout > 9;

wJvcount = 0;
wJcount = 0;
Jvcount = 0;
Jcount = 0;
wfcount = 0;
fcount = 0;

wJvundef = false;
wJundef = false;
Jvundef = false;
Jundef = false;

wJvapprox = false;
wJapprox = false;
Jvapprox = false;
Japprox = false;

iopts = opts;
iopts.FunValCheck = false;

mx = size(x, 1); % m is the number of individuals.
[mw, nobj] = size(w); % m is the number of individuals and nobj is the number of objectives.

if isempty(wJv) && isempty(wJ) && isempty(wJx) && isempty(Jv) && isempty(Jvx) && isempty(J) && isempty(Jx)
  wjacv();
else
  [wJvx, wJvcount, wJvundef, wJx, wJcount, wJundef, Jvx, Jvcount, Jvundef, Jx, Jcount, Jundef] = eval.wjvevalormult(...
    wJv, wJ, Jv, J, x, wJx, Jvx, Jx, w, v, iswarning, forcejac, opts);
  if wJvundef || wJundef || Jvundef || Jundef
    Jx = [];
    wjacv();
  end
end

  function [] = wjacv()
    wJvapprox = true;
    if forcejac
      if isempty(Jx)
        [Jx, ~, ~, Japprox, fx, fcount] = eval.jevalorfd([], f, x, fx, lb, ub, iswarning, opts);
      end
      [wJvx, wJvcount, ~, wJx, wJcount, ~, Jvx, Jvcount] = eval.wjvevalormult([], [], [], [], x, [], [], Jx, w, v, iswarning, false, opts);
    elseif mw > 1 && mx == 1
      [Jvx, ~, ~, Jvapprox, ~, ~, ~, ~, fx, fcount] = eval.jvevalormultorfd([], [], f, x, fx, [], v, lb, ub, iswarning, false, opts);
      wJvx = eval.wjevalormult([], [], x, Jvx, w, iswarning, false, opts);
    else
      topts = opts;
      topts.UseVectorized = true;

      if isempty(wfx)
        [wfx, wfcount, fx, fcount] = eval.wfeval(f, x, fx, w, false, opts);
      end
      [wJvx, ~, ~, Jvapprox, ~, ~, ~, ~, wfx, twfcount] = eval.jvevalormultorfd([], [], @(x) wfeval(x), x, wfx, [], v, lb, ub, iswarning, false, topts);
      wfcount = wfcount + twfcount;
      fcount = fcount + twfcount;
    end
  end

  function [wfy] = wfeval(y) 
    M = size(y, 1); % n * mx
    if mw == 1 || M == mx     
      wfy = eval.wfeval(f, y, [], w, false, iopts); % (M x 1)
    else
      rw = repmat(w, 1, M / mx)';
      rw = reshape(rw, nobj, M)';
      wfy = eval.wfeval(f, y, [], rw, false, iopts); % (M x 1)
    end
  end
end