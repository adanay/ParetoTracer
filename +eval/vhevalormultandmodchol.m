function [vHx, vHcount, vHundef,...
  Hx, Hcount, Hundef, Lx, Lmodif] = vhevalormultandmodchol(vH, H, x, Hx, v, iswarning, forcehess, opts)
% Vectorized Hessian multiply function evaluation plus a Modified Cholesky
% decomposition.
% vH = [v' * H1; v' * H2; ...; v' * Hnobj], where H1,..., Hnobj are modified
% to ensure they are positive definite.
% If vH is empty, it will perform the multiplication.
% If Lmodif = true, then Hx is not positive definite, thus Lx * Lx' needs
% to be computed.
% If Lmodif is empty, then the modification did not take place as e.g. the
% multiply function vH was provided or Hx was provided.
% Each row of x is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if
% opts.FunValCheck is true.
% forcehess will force one evaluation of the Hessian to be returned in Hx
% even if vH was passed and evaluated.
% The Cholesky modification will be performed only if the Hessian is
% defined and if it had to be computed. If the Hessian is passed as an
% argument it will be assumed to be positive definite and will not be
% modified.
% opts are the optimization options.
%
% The result is always of size (m x nobj x n) even if m = 1.
% Hx (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcehess = forcehess && nargout > 3;

vHundef = false;

Hcount = 0;
Hundef = false;

Lx = [];
Lmodif = [];

if isempty(vH)
  if isempty(Hx)
    [Hx, Hcount, Hundef, Lx, Lmodif] = eval.hevalandmodchol(H, x, iswarning, opts);
  end
  [vHx, vHcount] = eval.vhorvlltimes(Hx, Lx, Lmodif, v);
else
  [vHx, vHcount, vHundef] = eval.vheval(vH, x, v, iswarning, opts);
end

if forcehess && isempty(Hx)
  if nargout > 6
    [Hx, Hcount, Hundef, Lx, Lmodif] = eval.hevalandmodchol(H, x, iswarning, opts);
  else
    [Hx, Hcount, Hundef] = eval.hevalandmodchol(H, x, iswarning, opts);
  end
end
end