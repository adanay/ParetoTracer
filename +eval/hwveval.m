% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hwvx, Hwvcount, Hwvundef] = hwveval(Hwv, x, w, v, iswarning, opts)
% Vectorized Hessian multiply function evaluation.
% Hwv = (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x, w, and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x n) even if m = 1.

[mx, n] = size(x); % m is the number of individuals and n is the number of variables.
mw = size(w, 1); 
mv = size(v, 1); 
m = max([mx, mw, mv]);

Hwvcount = 0;
Hwvundef = false;

if m == 1 || opts.UseVectorized
  Hwvx = feval(Hwv, x, w, v);
  if m == 1
    Hwvx = Hwvx(:)';
  else
    Hwvx = Hwvx(:, :);
  end
  Hwvcount = m;
  if opts.FunValCheck 
    Hwvundef = val.checkval(Hwvx, 'Hwv', iswarning);
  end
else
  for i = 1 : m
    ix = min(i, mx);
    iw = min(i, mw);
    iv = min(i, mv);
    Hwvxi = feval(Hwv, x(ix, :), w(iw, :), v(iv, :));
    Hwvxi = Hwvxi(:)';
    
    if i == 1
      Hwvx = zeros(m, n);
    end
    
    Hwvx(i, :) = Hwvxi;
    Hwvcount = Hwvcount + 1;
    if opts.FunValCheck
      Hwvundef = Hwvundef || val.checkval(Hwvxi, 'Hwv', iswarning);
    end
  end
end
end