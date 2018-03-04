function [Hwx1, Hwcount, Hwundef, Hwapprox,... 
          Hx1, Hcount, Hundef, Happrox, Hident] =... 
  hwevalormultorqn(Hw, H, x1, Jx1, Hx1, x0, Jx0, Hx0, w, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation.
% Hw = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj.
% If Hw is empty, it will perform the multiplication.
% If H is also empty, a QN update will be performed.
% Each row of x1 (and x0) is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x1 and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and Hw(x1) or 
% H(x1) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be updated if the previous iteration values are provided.
% Otherwise it will be approximated to the identity matrix.
% + iswarning = false implies that an error will be thrown.
% forcehess will force one evaluation of the Hessian to be returned in Hx1 
% even if Hw was passed and evaluated.
% formident will force the formation of the identity to be returned in Hx1 
% in case there is no other choice than approximating Hx1 to the identity. 
% This may happen e.g. if the previous iteration values are not provided. 
% Otherwise Hx1 will be empty (but known to be the identity if Hident is true).
% opts are the optimization options.
%
% The result is always of size (m x n x n) even if m = 1.
% Jx1 (and Jx0) (if entered) is assumed to be of size (m x obj x n) even  
% if m = 1 (after being calculated using one of the vec eval functions).
% Hx1 (and Hx0) (if entered) is assumed to be of size (m x n x n x nobj)  
% even if m = 1 (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcehess = forcehess && nargout > 4;

Hwcount = 0;
Hcount = 0;

Hwundef = false;
Hwapprox = false;

Hundef = false;
Happrox = false;
Hident = false;

if isempty(Hw) && isempty(H) && isempty(Hx1)
  hessw();
else
  [Hwx1, Hwcount, Hwundef, Hx1, Hcount, Hundef] = eval.hwevalormult(Hw, H, x1, Hx1, w, iswarning, forcehess, opts);
  if Hwundef || Hundef
    if ~isempty(Jx1) && ~isempty(x0) && ~isempty(Jx0)
      hessw();
    else
      Hwapprox = true;
      Happrox = true;
      Hident = true;
      [Hwx1, ~, ~, ~, Hx1] = eval.hwevalormultorsd([], [], x1, [], w, iswarning, forcehess, formident, opts);
    end
  end
end

  function [] = hessw()
    Hwapprox = true;
    Happrox = true;
    [Hx1, ~, ~, ~, Hident] = eval.hevalorqn([], x1, Jx1, x0, Jx0, Hx0, iswarning, forcehess && formident, opts);
    if Hident
      Hwx1 = bsxfun(@times, sum(w, 2), utils.eyen([m n n], 2, 3)); % (m x n x n)
    else
      [Hwx1, tHwcount] = eval.hwevalormult([], [], x1, Hx1, w, iswarning, false, opts);
      Hwcount = Hwcount + tHwcount;
    end
  end
end
