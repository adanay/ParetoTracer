% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hwx, Hwcount, Hwundef, Lwx, Lwmodif] = hwevalandmodchol(Hw, x, w, iswarning, opts)
% Vectorized Hessian multiply function evaluation plus a modified Cholesky 
% decomposition.
% Hw = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj, where the resulting 
% weighted Hessian is modified to ensure it is positive definite.
% If Lwmodif = true, then Hwx is not positive definite, thus Lwx * Lwx'   
% needs to be computed.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true.
% The Cholesky modification will be performed only if the weighted Hessian 
% is defined.
% opts are the optimization options.
%
% The result is always of size (m x n x n) even if m = 1.

Lwx = [];
Lwmodif = [];

[Hwx, Hwcount, Hwundef] = eval.hweval(Hw, x, w, iswarning, opts);
if ~Hwundef && nargout > 3
  [Lwx, Lwmodif] = utils.modcholn(Hwx, 2, 3);
end
end