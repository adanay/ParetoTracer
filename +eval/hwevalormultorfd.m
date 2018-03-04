function [Hwx, Hwcount, Hwundef, Hwapprox,... 
          Hx, Hcount, Hundef, Happrox, Hident,...
          wfx, wfcount, fx, fcount,...
          wJx, wJcount, wJundef, Jx, Jcount, Jundef] =... 
  hwevalormultorfd(Hw, H, f, wJ, J, x, wfx, fx, wJx, Jx, Hx, w, lb, ub, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation.
% Hw = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj.
% If Hw is empty, it will perform the multiplication.
% If H is also empty, it will use finite differences (FD).
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and Hw(x) or 
% H(x) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be computed by FD if the Jacobian or the objective functions 
% are available. Otherwise it will be approximated to the identity.
% + iswarning = false implies that an error will be thrown.
% In case that FD is applied with the Jacobian, the value of the Jacobian 
% is the one that will be checked if opts.FunValCheck is true. If also
% iswarning = true, then FD will be re-computed using the objective 
% function f if available.
% forcehess will force one evaluation of the Hessian to be returned in Hx 
% even if Hw was passed and evaluated.
% formident will force the formation of the identity to be returned in Hx 
% in case there is no other choice than approximating Hx to the identity. 
% Otherwise Hx will be empty (but known to be the identity if Hident is true).
% opts are the optimization options.
%
% The result is always of size (m x n x n) even if m = 1.
% fx (if entered) is assumed to be of size (m x obj) even if m = 1 
% (after being calculated using one of the vec eval functions).
% wfx (if entered) is assumed to be of size (m x 1) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% wJx (if entered) is assumed to be of size (m x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Hx (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcehess = forcehess && nargout > 4;

Hwcount = 0;
Hwundef = false;
Hwapprox = false;

Hcount = 0;
Hundef = false;
Happrox = false;
Hident = false;

wfcount = 0;
fcount = 0;

wJcount = 0;
wJundef = false;
Jcount = 0;
Jundef = false;

[mx, n] = size(x); % m is the number of individuals and n is the number of variables.
[mw, nobj] = size(w); % nobj is the number of objectives.
m = max(mx, mw);

iopts = opts;
iopts.FunValCheck = false;

if isempty(Hw) && isempty(H) && isempty(Hx)
  hessw();
else
  [Hwx, Hwcount, Hwundef, Hx, Hcount, Hundef] = eval.hwevalormult(Hw, H, x, Hx, w, iswarning, forcehess, opts);
  if Hwundef || Hundef
    if ~isempty(f) || ~isempty(J) || ~isempty(wJ) && ~(forcehess || mw > 1 && mx == 1)
      Hx = [];
      hessw();
    else
      Hwapprox;
      Happrox = true;
      Hident = true;
      [Hwx, ~, ~, ~, Hx] = eval.hwevalormultorsd([], [], x, [], w, iswarning, forcehess, formident, opts);
    end
  end
end

  function [] = hessw()
    Hwapprox = true;
    if forcehess || mw > 1 && mx == 1
      if isempty(Hx)
        [Hx, ~, ~, Happrox, Hident, fx, fcount, Jx, Jcount, Jundef] = eval.hevalorfd([], f, J, x, fx, Jx, lb, ub, iswarning, forcehess && formident, opts);
      end
      if Hident
        Hwx = bsxfun(@times, sum(w, 2), utils.eyen([m n n], 2, 3)); % (m x n x n)
      else
        [Hwx, tHwcount] = eval.hwevalormult([], [], x, Hx, w, iswarning, false, opts);
        Hwcount = Hwcount + tHwcount;
      end
    else
      if isempty(J) && isempty(wJ)
        hessw_f();
      else
        hessw_j();
      end
    end
  end

  function [] = hessw_j()
    topts = opts;
    topts.UseVectorized = true;
    
    if isempty(wJx)
      [wJx, wJcount, wJundef, Jx, Jcount, Jundef] = eval.wjevalormult(wJ, J, x, Jx, w, true, false, opts);
    end
    [Hwx, ~, ~, ~, ~, ~, ~, wJx, twJcount, twJundef] = eval.hevalorfd([], [], @(x) wjeval(x), x, [], permute(wJx, [1 3 2]), lb, ub, iswarning, true, topts);
    wJx = permute(wJx, [1 3 2]);
    wJcount = wJcount + twJcount;
    wJundef = wJundef || twJundef;
    if isempty(wJ)
      Jcount = Jcount + twJcount;
      Jundef = Jundef || twJundef;
    end
    if Jundef || wJundef
      if ~isempty(f)
        hessw_f();
      else
        Hwx = eval.hwevalormultorsd([], [], x, [], w, iswarning, false, false, opts);
      end
    end
  end

  function [] = hessw_f()
    topts = opts;
    topts.UseVectorized = true;
    
    if isempty(wfx)
      [wfx, wfcount, fx, fcount] = eval.wfeval(f, x, fx, w, false, opts);
    end
    [Hwx, ~, ~, ~, ~, wfx, twfcount] = eval.hevalorfd([], @(x) wfeval(x), [], x, wfx, [], lb, ub, iswarning, true, topts);
    wfcount = wfcount + twfcount;
    fcount = fcount + twfcount;
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

  function [wJy] = wjeval(y)
    M = size(y, 1); % n * mx
    if mw == 1 || M == mx
      wJy = eval.wjevalormult(wJ, J, y, [], w, false, false, iopts); % (M x n)
    else
      rw = repmat(w, 1, M / mx)';
      rw = reshape(rw, nobj, M)';
      wJy = eval.wjevalormult(wJ, J, y, [], rw, false, false, iopts); % (M x n)
    end   
    wJy = permute(wJy, [1 3 2]); % (M x 1 x n)
  end
end