% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Jx, Jcount, Jundef] = jeval(J, x, iswarning, opts)
% Vectorized Jacobian function evaluation.
% Each row of x is considered an individual to be evaluated.
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x nobj x n) even if m = 1.
% The result of J(x) is assumed to be (nobj x n) if m = 1.

Jcount = 0;
Jundef = false;

[m, n] = size(x); % m is the number of individuals and n is the number of variables.

if m == 1 || opts.UseVectorized
  Jx = feval(J, x);
  if m == 1
    Jx = Jx(:, :);
    Jx = shiftdim(Jx, -1); % (1 x nobj x n)
  else
    Jx = Jx(:, :, :);
  end
  Jcount = m;
  if opts.FunValCheck 
    Jundef = val.checkval(Jx, 'J', iswarning);
  end
else
  for i = 1 : m
    Jxi = feval(J, x(i, :));
    Jxi = Jxi(:, :);
    
    if i == 1
      nobj = size(Jxi, 1);
      Jx = zeros(m, nobj, n);
    end
    
    Jx(i, :, :) = Jxi;
    Jcount = Jcount + 1;
    if opts.FunValCheck
      Jundef = Jundef || val.checkval(Jxi, 'J', iswarning);
    end
  end
end
end