function [Hwvx1, Hwvcount, Hwvundef, Hwvapprox,... 
          Hwx1, Hwcount, Hwundef,... 
          vHx1, vHcount, vHundef,... 
          Hx1, Hcount, Hundef, Happrox, Hident] =... 
  hwvevalormultorqn(Hwv, Hw, vH, H, x1, Jx1, Hwx1, vHx1, Hx1, x0, Jx0, Hx0, w, v, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation.
% Hwv = (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v.
% If Hwv is empty, it will perform the multiplication.
% If H is also empty, a QN update will be performed.
% Each row of x1 (and x0) is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x1, w, and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and Hwv(x1) or 
% Hw(x1), or vH(1x), or H(1x) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be updated if the previous iteration values are provided.
% Otherwise it will be approximated to the identity matrix.
% + iswarning = false implies that an error will be thrown.
% forcehess will force one evaluation of the Hessian to be returned in Hx 
% even if Hwv was passed and evaluated.
% formident will force the formation of the identity to be returned in Hx1 
% in case there is no other choice than approximating Hx1 to the identity. 
% This may happen e.g. if the previous iteration values are not provided. 
% Otherwise Hx1 will be empty (but known to be the identity if Hident is true).
% opts are the optimization options.
%
% The result is always of size (m x n) even if m = 1.
% Jx1 (and Jx0) (if entered) is assumed to be of size (m x obj x n) even  
% if m = 1 (after being calculated using one of the vec eval functions).
% Hx1 (and Hx0) (if entered) is assumed to be of size (m x n x n x nobj)  
% even if m = 1 (after being calculated using one of the vec eval functions).
% Hwx1 (if entered) is assumed to be of size (m x n x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% vHx1 (if entered) is assumed to be of size (m x nobj x n) even if m = 1 
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
Hident = false;

if isempty(Hwv) && isempty(Hw) && isempty(Hwx1) && isempty(vH) && isempty(vHx1) && isempty(H) && isempty(Hx1)
  hesswv();
else
  [Hwvx1, Hwvcount, Hwvundef, Hwx1, Hwcount, Hwundef, vHx1, vHcount, vHundef, Hx1, Hcount, Hundef] = eval.hwvevalormult(...
    Hwv, Hw, vH, H, x1, Hwx1, vHx1, Hx1, w, v, iswarning, forcehess, opts);
  if Hwvundef || Hwundef || vHundef || Hundef
    if ~isempty(Jx1) && ~isempty(x0) && ~isempty(Jx0)
      hesswv();
    else
      Hwvapprox = true;
      Happrox = true;
      Hident = true;
      [Hwvx1, ~, ~, ~, Hwx1, ~, ~, vHx1, ~, ~, Hx1] = eval.hwvevalormultorsd([], [], [], [], x1, [], [], [], w, v, iswarning, forcehess, formident, opts);
    end
  end
end

  function [] = hesswv()
    Hwvapprox = true;
    Happrox = true;
    
    [Hx1, ~, ~, ~, Hident] = eval.hevalorqn([], x1, Jx1, x0, Jx0, Hx0, iswarning, forcehess && formident, opts);
    
    if Hident
      Hwvx1 = bsxfun(@times, sum(w, 2), v); % (m x n)
      if size(Hwvx1, 1) < m
        Hwvx1 = repmat(Hwvx1, m, 1);
      end
    else
      [Hwvx1, tHwvcount, ~, Hwx1, tHwcount, ~, vHx1, tvHcount] = eval.hwvevalormult([], [], [], [], x1, [], [], Hx1, w, v, iswarning, false, opts);
      Hwvcount = Hwvcount + tHwvcount;
      Hwcount = Hwcount + tHwcount;
      vHcount = vHcount + tvHcount;
    end
  end
end