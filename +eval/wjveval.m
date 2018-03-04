function [wJvx, wJvcount, wJvundef] = wjveval(wJv, x, w, v, iswarning, opts)
% Vectorized Jacobian multiply function evaluation.
% wJv = w' * J * v.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x, w, and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x 1).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

wJvcount = 0;
wJvundef = false;

mx = size(x, 1); % m is the number of individuals.
mw = size(w, 1); 
mv = size(v, 1); 
m = max([mx, mw, mv]);

if m == 1 || opts.UseVectorized
  wJvx = feval(wJv, x, w, v);
  wJvx = wJvx(:, :);
  wJvcount = m;
  if opts.FunValCheck 
    wJvundef = val.checkval(wJvx, 'wJv', iswarning);
  end
else
  for i = 1 : m
    ix = min(i, mx);
    iw = min(i, mw);
    iv = min(i, mv);
    wJvxi = feval(wJv, x(ix, :), w(iw, :), v(iv, :));
    wJvxi = wJvxi(:)';
    
    if i == 1
      wJvx = zeros(m, 1);
    end
    
    wJvx(i) = wJvxi;
    wJvcount = wJvcount + 1;
    if opts.FunValCheck
      wJvundef = wJvundef || val.checkval(wJvxi, 'wJv', iswarning);
    end
  end
end
end