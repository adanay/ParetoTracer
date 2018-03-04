function [Hwx, Hwcount, Hwundef] = hweval(Hw, x, w, iswarning, opts)
% Vectorized Hessian multiply function evaluation.
% Hw = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x n x n) even if m = 1.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

[mx, n] = size(x); % m is the number of individuals and n is the number of variables.
mw = size(w, 1); 
m = max(mx, mw);

Hwcount = 0;
Hwundef = false;

if m == 1 || opts.UseVectorized
  Hwx = feval(Hw, x, w);
  if m == 1
    Hwx = Hwx(:, :);
    Hwx = shiftdim(Hwx, -1); % (1 x n x n)
  else
    Hwx = Hwx(:, :, :);
  end
  Hwcount = m;
  if opts.FunValCheck 
    Hwundef = val.checkval(Hwx, 'Hw', iswarning);
  end
else
  for i = 1 : m
    ix = min(i, mx);
    iw = min(i, mw);
    Hwxi = feval(Hw, x(ix, :), w(iw, :));
    Hwxi = Hwxi(:, :);
    
    if i == 1
      Hwx = zeros(m, n, n);
    end
    
    Hwx(i, :, :) = Hwxi;
    Hwcount = Hwcount + 1;
    if opts.FunValCheck
      Hwundef = Hwundef || val.checkval(Hwxi, 'Hw', iswarning);
    end
  end
end
end