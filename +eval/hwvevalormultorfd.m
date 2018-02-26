% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hwvx, Hwvcount, Hwvundef, Hwvapprox,... 
          Hwx, Hwcount, Hwundef,... 
          vHx, vHcount, vHundef, vHapprox,... 
          Hx, Hcount, Hundef, Happrox, Hident,...
          wfx, wfcount, fx, fcount,...
          wJx, wJcount, wJundef, Jx, Jcount, Jundef] =... 
  hwvevalormultorfd(Hwv, Hw, vH, H, f, wJ, J, x, wfx, fx, wJx, Jx, Hwx, vHx, Hx, w, v, lb, ub, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation.
% Hwv = (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v.
% If Hwv is empty, it will perform the multiplication.
% If H is also empty, it will use finite differences (FD).
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x, w, and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and Hwv(x) or 
% Hw(x), or vH(x), or H(x) is undef, then
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
% The result is always of size (m x n) even if m = 1.
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
% Hwx (if entered) is assumed to be of size (m x n x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% vHx (if entered) is assumed to be of size (m x nobj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).

forcehess = forcehess && nargout > 11;

Hwvcount = 0;
Hwvundef = 0;
Hwvapprox = false;

Hwcount = 0;
Hwundef = 0;

vHcount = 0;
vHundef = 0;
vHapprox = false;

Hcount = 0;
Hundef = 0;
Happrox = false;
Hident = false;

wfcount = 0;
fcount = 0;

wJcount = 0;
wJundef = 0;

Jcount = 0;
Jundef = 0;

mx = size(x, 1); % m is the number of individuals and n is the number of variables.
[mw, nobj] = size(w); % nobj is the number of objectives.
mv = size(v, 1);
m = max([mx, mw, mv]);

iopts = opts;
iopts.FunValCheck = false;

if isempty(Hwv) && isempty(Hw) && isempty(Hwx) && isempty(vH) && isempty(vHx) && isempty(H) && isempty(Hx)
  hesswv();
else
  [Hwvx, Hwvcount, Hwvundef, Hwx, Hwcount, Hwundef, vHx, vHcount, vHundef, Hx, Hcount, Hundef] = eval.hwvevalormult(...
    Hwv, Hw, vH, H, x, Hwx, vHx, Hx, w, v, iswarning, forcehess, opts);
  if Hwvundef || Hwundef || vHundef || Hundef
    if ~isempty(f) || ~isempty(J) || ~isempty(wJ) && ~(forcehess || mw > 1 && mx == 1)
      Hx = [];
      hesswv();
    else
      Hwvapprox;
      Happrox = true;
      Hident = true;
      [Hwvx, ~, ~, ~, Hwx, ~, ~, vHx, ~, ~, Hx] = eval.hwvevalormultorsd([], [], [], [], x, [], [], [], w, v, iswarning, forcehess, formident, opts);
    end
  end
end

  function [] = hesswv()
    Hwvapprox = true;
    if forcehess
      if isempty(Hx)
        [Hx, ~, ~, Happrox, Hident, fx, fcount, Jx, Jcount, Jundef] = eval.hevalorfd([], f, J, x, fx, Jx, lb, ub, iswarning, forcehess && formident, opts);
      end
      if Hident
        Hwvx = bsxfun(@times, sum(w, 2), v); % (m x n)
        if size(Hwvx, 1) < m
          Hwvx = repmat(Hwvx, m, 1);
        end
      else
        [Hwvx, tHwvcount, ~, Hwx, tHwcount, ~, vHx, tvHcount] = eval.hwvevalormult([], [], [], [], x, [], [], Hx, w, v, iswarning, false, opts);
        Hwvcount = Hwvcount + tHwvcount;
        Hwcount = Hwcount + tHwcount;
        vHcount = vHcount + tvHcount;
      end
    elseif mw > 1 && mx == 1
      [vHx, ~, ~, vHapprox, ~, ~, ~, ~, ~, fx, fcount, Jx, Jcount, Jundef] = eval.vhevalormultorfd([], [], f, J, x, fx, Jx, [], v, lb, ub, iswarning, false, false, opts); % (m x nobj x n)
      Hwvx = eval.wjevalormult([], [], x, vHx, w, iswarning, false, opts); % wjevalormult is ok
    else
      if isempty(J) && isempty(wJ)
        hesswv_f();
      else
        hesswv_j();
      end
      Hwvx = permute(Hwvx, [1 3 2]);
    end
  end

  function [] = hesswv_j()
    topts = opts;
    topts.UseVectorized = true;
    
    if isempty(wJx)
      [wJx, wJcount, wJundef, Jx, Jcount, Jundef] = eval.wjevalormult(wJ, J, x, Jx, w, true, false, opts);
    end
    [Hwvx, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, wJx, twJcount, twJundef] = eval.vhevalormultorfd([], [], [], @(x) wjeval(x), x, [], permute(wJx, [1 3 2]), [], v, lb, ub, iswarning, false, false, topts);
    wJx = permute(wJx, [1 3 2]);
    wJcount = wJcount + twJcount;
    wJundef = wJundef || twJundef;
    if isempty(wJ)
      Jcount = Jcount + twJcount;
      Jundef = Jundef || twJundef;
    end
    if Jundef || wJundef
      if ~isempty(f)
        hesswv_f();
      else
        Hwvx = eval.hwvevalormultorsd([], [], [], [], x, [], [], [], w, v, iswarning, false, false, opts);
      end
    end
  end

  function [] = hesswv_f()
    topts = opts;
    topts.UseVectorized = true;
    
    if isempty(wfx)
      [wfx, wfcount, fx, fcount] = eval.wfeval(f, x, fx, w, false, opts);
    end
    [Hwvx, ~, ~, ~, ~, ~, ~, ~, ~, wfx, twfcount] = eval.vhevalormultorfd([], [], @(x) wfeval(x), [], x, wfx, [], [], v, lb, ub, iswarning, false, true, topts);
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