function [vHx, vHcount, vHundef] = vheval(vH, x, v, iswarning, opts)
% Vectorized Hessian multiply function evaluation.
% vH = [v' * H1; v' * H2; ...; v' * Hnobj].
% Each row of x is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x nobj x n) even if m = 1.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

[mx, n] = size(x); % m is the number of individuals and n is the number of variables.
mv = size(v, 1);
m = max(mx, mv);

vHcount = 0;
vHundef = false; 

if m == 1 || opts.UseVectorized
  vHx = feval(vH, x, v);
  if m == 1
    vHx = vHx(:, :);
    vHx = shiftdim(vHx, -1); % (1 x nobj x n)
  else
    vHx = vHx(:, :, :);
  end
  vHcount = m;
  if opts.FunValCheck 
    vHundef = val.checkval(vHx, 'vH', iswarning);
  end
else
  for i = 1 : m
    ix = min(i, mx);
    iv = min(i, mv);
    vHxi = feval(vH, x(ix, :), v(iv, :));
    vHxi = vHxi(:, :);
    
    if i == 1
      nobj = size(vHxi, 1);
      vHx = zeros(m, nobj, n);
    end
    
    vHx(i, :, :) = vHxi;
    vHcount = vHcount + 1;
    if opts.FunValCheck
      vHundef = vHundef || val.checkval(vHxi, 'vH', iswarning);
    end
  end
end
end