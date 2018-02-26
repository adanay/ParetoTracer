% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [wJx, wJcount, wJundef] = wjeval(wJ, x, w, iswarning, opts)
% Vectorized Jacobian multiply function evaluation.
% wJ = w' * J.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x n) even if m = 1.

wJcount = 0;
wJundef = false;

[mx, n] = size(x); % m is the number of individuals and n is the number of variables.
mw = size(w, 1); 
m = max(mx, mw);

if m == 1 || opts.UseVectorized
  wJx = feval(wJ, x, w);
  if m == 1
    wJx = wJx(:)'; % (1 x n)
  else
    wJx = wJx(:, :);
  end
  wJcount = m;
  if opts.FunValCheck 
    wJundef = val.checkval(wJx, 'wJ', iswarning);
  end
else
  for i = 1 : m
    ix = min(i, mx);
    iw = min(i, mw);
    wJxi = feval(wJ, x(ix, :), w(iw, :));
    wJxi = wJxi(:)';
    
    if i == 1
      wJx = zeros(m, n);
    end
    
    wJx(i, :) = wJxi;
    wJcount = wJcount + 1;
    if opts.FunValCheck
      wJundef = wJundef || val.checkval(wJxi, 'wJ', iswarning);
    end
  end
end
end