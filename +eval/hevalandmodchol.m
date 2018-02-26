% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hx, Hcount, Hundef, Lx, Lmodif] = hevalandmodchol(H, x, iswarning, opts)
% Vectorized Hessian function evaluation plus a Modified Cholesky 
% decomposition.
% If Lmodif = true, then Hx is not positive definite, thus Lx * Lx' needs  
% to be computed.
% Each row of x is considered an individual to be evaluated.
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true.
% The Cholesky modification will be performed only if the Hessian is
% defined.
% opts are the optimization options.
%
% The result is always of size (m x n x n x nobj) even if m = 1.
% The result of H(x) is assumed to be (n x n x nobj) if m = 1.

Lx = [];
Lmodif = [];

[Hx, Hcount, Hundef] = eval.heval(H, x, iswarning, opts);
if ~Hundef && nargout > 3
  [Lx, Lmodif] = utils.modcholn(Hx, 2, 3);
end
end