% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hwx, Hwcount, Hwundef, Hwapprox, Lwx, Lwmodif,... 
          Hx, Hcount, Hundef, Happrox, Hident,...
          wfx, wfcount, fx, fcount,...
          wJx, wJcount, wJundef, Jx, Jcount, Jundef] =... 
  hwevalormultorfdandmodchol(Hw, H, f, wJ, J, x, wfx, fx, wJx, Jx, Hx, w, lb, ub, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation plus a Modified Cholesky
% decomposition.
% Hw = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj, where the resulting 
% weighted Hessian is modified to ensure it is positive definite.
% Note that the individual Hessians H1,..., Hnobj will never be modified. 
% If Hw is empty, it will perform the multiplication.
% If H is also empty, it will use finite differences (FD).
% A Modified Cholesky decomposition is performed afterwards.
% If Lwmodif = true, then Hwx is not positive definite, thus Lwx * Lwx'   
% needs to be computed.
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

Lwx = [];
Lwmodif = [];

[Hwx, Hwcount, Hwundef, Hwapprox,... 
 Hx, Hcount, Hundef, Happrox, Hident,...
 wfx, wfcount, fx, fcount,...
 wJx, wJcount, wJundef, Jx, Jcount, Jundef] = eval.hwevalormultorfd(Hw, H, f, wJ, J, x, wfx, fx, wJx, Jx, Hx, w, lb, ub, iswarning, forcehess, formident, opts);

if nargout > 4
  [Lwx, Lwmodif] = utils.modcholn(Hwx, 2, 3);
end

if nargout > 4
  if Hident
    tol = eps;
    s = sum(w, 2);
    if s > tol
      Lwmodif = false;
      Lwx = bsxfun(@times, sqrt(s), utils.eyen([m n n], 2, 3)); 
    else
      [Lwx, Lwmodif] = utils.modcholn(Hwx, 2, 3);
    end
  else
    [Lwx, Lwmodif] = utils.modcholn(Hwx, 2, 3);
  end
end
end