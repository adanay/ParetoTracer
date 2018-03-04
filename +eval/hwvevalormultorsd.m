function [Hwvx, Hwvcount, Hwvundef, Hwvapprox,... 
          Hwx, Hwcount, Hwundef,...
          vHx, vHcount, vHundef,...
          Hx, Hcount, Hundef, Happrox] =... 
  hwvevalormultorsd(Hwv, Hw, vH, H, x, Hwx, vHx, Hx, w, v, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation.
% Hwv = (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v.
% If Hwv is empty, it will perform the multiplication.
% If H is also empty, the Hessian will be approximated to the identity matrix.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x, w, and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and Hwv(x) or 
% Hw(x), or vH(x), or H(x) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be approximated to the identity matrix.
% + iswarning = false implies that an error will be thrown.
% forcehess will force one evaluation of the Hessian to be returned in Hx 
% even if Hw was passed and evaluated.
% formident will force the formation of the identity to be returned in Hx 
% if it is approximated. Otherwise Hx will be empty (but known to be the 
% identity if Happrox is true).
% opts are the optimization options.
%
% The result is always of size (m x n) even if m = 1.
% Hx (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Hwx (if entered) is assumed to be of size (m x n x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% vHx (if entered) is assumed to be of size (m x nobj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcehess = forcehess && nargout > 10;

Hwvcount = 0;
Hwvundef = false;
Hwvapprox = false;

Hwcount = 0;
Hwundef = false;

vHcount = 0;
vHundef = false;

Hcount = 0;
Hundef = false;
Happrox = false;

[mx, n] = size(x); % m is the number of individuals and n is the number of variables.
[mw, nobj] = size(w); % nobj is the number of objectives.
mv = size(v, 1);
m = max([mx, mw, mv]);

if isempty(Hwv) && isempty(Hw) && isempty(Hwx) && isempty(vH) && isempty(vHx) && isempty(H) && isempty(Hx)
  hesswv();
else
  [Hwvx, Hwvcount, Hwvundef, Hwx, Hwcount, Hwundef, vHx, vHcount, vHundef, Hx, Hcount, Hundef] = eval.hwvevalormult(...
    Hwv, Hw, vH, H, x, Hwx, vHx, Hx, w, v, iswarning, forcehess, opts);
  if Hwvundef || Hwundef || vHundef || Hundef
    Hx = [];
    hesswv();
  end
end

  function [] = hesswv()
    Hwvapprox = true;
    Happrox = true;
    
    Hwvx = bsxfun(@times, sum(w, 2), v); % (m x n)
    if size(Hwvx, 1) < m
      Hwvx = repmat(Hwvx, m, 1);
    end
    
    if forcehess && isempty(Hx)
      if formident
        Hx = utils.eyen([mx n n nobj], 2, 3);
      end
    end
  end
end