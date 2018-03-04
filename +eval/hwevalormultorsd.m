function [Hwx, Hwcount, Hwundef, Hwapprox,... 
          Hx, Hcount, Hundef, Happrox] =...
  hwevalormultorsd(Hw, H, x, Hx, w, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation.
% Hw = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj.
% If Hw is empty, it will perform the multiplication.
% If H is also empty, the Hessian will be approximated to the identity matrix.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and Hw(x) or 
% H(x) is undef, then
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
% The result is always of size (m x n x n) even if m = 1.
% Hx (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).

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

[mx, n] = size(x); % m is the number of individuals and n is the number of variables.
[mw, nobj] = size(w); % nobj is the number of objectives.
m = max(mx, mw);

if isempty(Hw) && isempty(H) && isempty(Hx)
  hessw();
else
  [Hwx, Hwcount, Hwundef, Hx, Hcount, Hundef] = eval.hwevalormult(Hw, H, x, Hx, w, iswarning, forcehess, opts);
  if Hwundef || Hundef
    Hx = [];
    hessw();
  end
end

  function [] = hessw()
    Hwapprox = true; 
    Happrox = true;
    
    Hwx = bsxfun(@times, sum(w, 2), utils.eyen([m n n], 2, 3)); % (m x n x n)
    
    if forcehess && isempty(Hx)
      if formident
        Hx = utils.eyen([m n n nobj], 2, 3);
      end
    end
  end
end
