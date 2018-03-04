function [Hx, Hcount, Hundef] = heval(H, x, iswarning, opts)
% Vectorized Hessian function evaluation.
% Each row of x is considered an individual to be evaluated.
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x n x n x nobj) even if m = 1.
% The result of H(x) is assumed to be (n x n x nobj) if m = 1.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

Hcount = 0;
Hundef = false;

[m, n] = size(x); % m is the number of individuals and n is the number of variables.

if m == 1 || opts.UseVectorized
  Hx = feval(H, x);
  if m == 1
    Hx = Hx(:, :, :);
    Hx = shiftdim(Hx, -1); % (1 x n x n x nobj)
  else
    Hx = Hx(:, :, :, :);
  end
  Hcount = m;
  if opts.FunValCheck 
    Hundef = val.checkval(Hx, 'H', iswarning);
  end
else
  for i = 1 : m
    Hxi = feval(H, x(i, :));
    Hxi = Hxi(:, :, :);
    
    if i == 1
      nobj = size(Hxi, 3);
      Hx = zeros(m, n, n, nobj);
    end
    
    Hx(i, :, :, :) = Hxi;
    Hcount = Hcount + 1;
    if opts.FunValCheck
      Hundef = Hundef || val.checkval(Hxi, 'H', iswarning);
    end
  end
end
end