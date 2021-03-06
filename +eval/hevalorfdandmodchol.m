function [Hx, Hcount, Hundef, Happrox, Hident, Lx, Lmodif,...
          fx, fcount, Jx, Jcount, Jundef, nobj] = hevalorfdandmodchol(H, f, J, x, fx, Jx, lb, ub, iswarning, formident, opts)
% Vectorized Hessian function evaluation plus a modified Cholesky 
% decomposition.
% If H is empty, finite differences (FD) are used.
% A Modified Cholesky decomposition is performed afterwards.
% If Lmodif = true, then Hx is not positive definite, thus Lx * Lx' needs  
% to be computed.
% Each row of x is considered an individual to be evaluated.
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and H(x) is 
% undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be computed by FD if the Jacobian or the objective functions 
% are available. Otherwise the Hessian will be approximated to the identity 
% matrix.
% + iswarning = false implies that an error will be thrown.
% In case that FD is applied with the Jacobian, the value of the Jacobian 
% is the one that will be checked if opts.FunValCheck is true. If also
% iswarning = true, then FD will be re-computed using the objective 
% function f if available.
% formident will force the formation of the identity to be returned in Hx 
% (and Lx) in case there is no other choice than approximating Hx to the 
% identity. Otherwise Hx (and Lx) will be empty (but known to be the 
% identity if Hident is true).
% opts are the optimization options.
%
% The result is always of size (m x n x n x nobj) even if m = 1.
% The result of H(x) is assumed to be (n x n x nobj) if m = 1.
% fx (if entered) is assumed to be of size (m x obj) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

Lx = [];
Lmodif = [];

[Hx, Hcount, Hundef, Happrox, Hident,...
 fx, fcount, Jx, Jcount, Jundef, nobj] = eval.hevalorfd(H, f, J, x, fx, Jx, lb, ub, iswarning, formident, opts);

if nargout > 5
  if Hident
    Lmodif = false;
    if formident
      Lx = Hx;
    end
  else
    [Lx, Lmodif] = utils.modcholn(Hx, 2, 3);
  end
end
end