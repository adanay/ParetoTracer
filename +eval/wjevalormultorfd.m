function [wJx, wJcount, wJundef, wJapprox,...
         Jx, Jcount, Jundef, Japprox,...
         wfx, wfcount, fx, fcount] = wjevalormultorfd(...
  wJ, J, f, x, wfx, fx, Jx, w, lb, ub, iswarning, forcejac, opts)
% Vectorized Jacobian multiply function evaluation.
% wJ = w' * J.
% If wJ is empty, it will perform the multiplication.
% If J is also empty, it will use finite differences (FD).
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and wJ(x) or 
% J(x) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Jacobian will be computed by FD (the objective function must be available).
% + iswarning = false implies that an error will be thrown.
% forcejac will force one evaluation of the Jacobian to be returned in Jx 
% even if wJ was passed and evaluated.
% opts are the optimization options.
%
% The result is always of size (m x n) even if m = 1.
% fx (if entered) is assumed to be of size (m x obj) even if m = 1 
% (after being calculated using one of the vec eval functions).
% wfx (if entered) is assumed to be of size (m x 1) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcejac = forcejac && nargout > 3;

wJcount = 0;
Jcount = 0;
fcount = 0;
wfcount = 0;

wJundef = false;
wJapprox = false;

Jundef = false;
Japprox = false;

iopts = opts;
iopts.FunValCheck = false;

mx = size(x, 1); % m is the number of individuals.
[mw, nobj] = size(w); % nobj is the number of objectives.

if isempty(wJ) && isempty(J) && isempty(Jx)
  wjac();
else
  [wJx, wJcount, wJundef, Jx, Jcount, Jundef] = eval.wjevalormult(wJ, J, x, Jx, w, iswarning, forcejac, opts);
  if wJundef || Jundef
    Jx = [];
    wjac();
  end
end

  function [] = wjac()
    wJapprox = true;
    if forcejac || (mw > 1 && mx == 1)
      if isempty(Jx)
        [Jx, ~, ~, Japprox, fx, fcount] = eval.jevalorfd([], f, x, fx, lb, ub, iswarning, opts);
      end
      [wJx, wJcount] = eval.wjevalormult([], [], x, Jx, w, iswarning, false, opts);
    else
      topts = opts;
      topts.UseVectorized = true;
      
      if isempty(wfx)
        [wfx, wfcount, fx, fcount] = eval.wfeval(f, x, fx, w, false, opts);
      end
      [wJx, ~, ~, ~, wfx, twfcount] = eval.jevalorfd([], @(x) wfeval(x), x, wfx, lb, ub, iswarning, topts); % (m x 1 x n)
      wJx = permute(wJx, [1 3 2]); % (m x n)
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