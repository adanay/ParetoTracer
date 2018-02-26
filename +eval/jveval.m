% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Jvx, Jvcount, Jvundef] = jveval(Jv, x, v, iswarning, opts)
% Vectorized Jacobian multiply function evaluation.
% Jv = J * v.
% Each row of x is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x nobj) even if m = 1.
% The result of Jv(x) is assumed to be (nobj x 1) if m = 1.

Jvcount = 0;
Jvundef = false;

mx = size(x, 1); % m is the number of individuals.
mv = size(v, 1); 
m = max(mx, mv);

if m == 1 || opts.UseVectorized
  Jvx = feval(Jv, x, v);
  if m == 1
    Jvx = Jvx(:)'; % (1 x nobj)
  else
    Jvx = Jvx(:, :);
  end
  Jvcount = m;
  if opts.FunValCheck 
    Jvundef = val.checkval(Jvx, 'Jv', iswarning);
  end
else
  for i = 1 : m
    ix = min(i, mx);
    iv = min(i, mv);
    Jvxi = feval(Jv, x(ix, :), v(iv, :));
    Jvxi = Jvxi(:)';
    
    if i == 1
      nobj = length(Jvxi);
      Jvx = zeros(m, nobj);
    end
    
    Jvx(i, :) = Jvxi;
    Jvcount = Jvcount + 1;
    if opts.FunValCheck
      Jvundef = Jvundef || val.checkval(Jvxi, 'Jv', iswarning);
    end
  end
end
end