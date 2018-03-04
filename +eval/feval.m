function [fx, fcount, fundef] = feval(f, x, iswarning, opts)
% Vectorized objective function evaluation.
% Each row of x is considered an individual to be evaluated.
% iswarning defines whether to display a warning or an error if the
% function value is undefined. The function value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x nobj) even if m = 1.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

fcount = 0;
fundef = false;

m = size(x, 1); % no of individuals

if m == 1 || opts.UseVectorized
  fx = feval(f, x);
  if m == 1
    fx = fx(:)';
  else
    fx = fx(:, :);
  end
  fcount = m; 
  if opts.FunValCheck 
    fundef = val.checkval(fx, 'f', iswarning);
  end
else
  for i = 1 : m
    fxi = feval(f, x(i, :));
    fxi = fxi(:)';
    
    if i == 1
      nobj = length(fxi);
      fx = zeros(m, nobj);
    end
    
    fx(i, :) = fxi;
    fcount = fcount + 1; 
    if opts.FunValCheck
      fundef = fundef || val.checkval(fxi, 'f', iswarning);
    end
  end
end
end

