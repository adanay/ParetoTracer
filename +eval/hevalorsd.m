% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hx, Hcount, Hundef, Happrox] = hevalorsd(H, x, nobj, iswarning, formident, opts)
% Vectorized Hessian function evaluation.
% If H is empty, the Hessian will be approximated to the identity matrix.
% Each row of x is considered an individual to be evaluated.
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and H(x) is 
% undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be approximated to the identity matrix.
% + iswarning = false implies that an error will be thrown.
% formident will force the formation of the identity to be returned in Hx 
% if it is approximated. Otherwise Hx will be empty (but known to be the 
% identity if Happrox is true).
% opts are the optimization options.
%
% The result is always of size (m x n x n x nobj) even if m = 1.
% The result of H(x) is assumed to be (n x n x nobj) if m = 1.

Hcount = 0;
Hundef = false;
Happrox = false;

[m, n] = size(x); % m is the number of individuals and n is the number of variables.

if isempty(H)                
  hess();
else
  [Hx, Hcount, Hundef] = eval.heval(H, x, iswarning, opts);
  if Hundef
    hess();
  end
end

  function [] = hess()
    Happrox = true;
    Hx = [];
    if formident
      Hx = utils.eyen([m n n nobj], 2, 3);
    end
  end
end