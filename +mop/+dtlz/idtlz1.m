% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = idtlz1(sizes, doval)
% DTLZ1 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% Suggested nobj = 3, n = nobj + 4 = 7.

% problem sizes
[n, nobj] = val.valmopsizes(sizes, true);
if n < 2
  n = 2;
end
if nobj < 2
  nobj = 2;
end
if nobj > n
  nobj = n;
end

sizes = struct(...
  'variables', n,...
  'objectives', nobj);

if nargin < 2
  doval = false;
end

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset] = mop.dtlz.dtlz(sizes, @idtlz1f);
pfront = @idtlz1pfront;
opts.MOPName = 'IDTLZ1';

  function [fx] = idtlz1f(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    q = n - nobj + 1;
    
    y = x(:, 1 : nobj - 1);
    z = x(:, nobj : n);
    
    gz = 100 * (q + sum((z - 0.5).^2 - cos(20 * pi * (z - 0.5)), 2));  
    
    hy = cumprod([ones(m, 1), y], 2);
    hy = fliplr(hy);
    
    ry = [ones(m, 1), 1 - y(:, end : -1 : 1)];
    
    fx = 0.5 * bsxfun(@minus, 1 + gz, bsxfun(@times, 1 + gz, hy .* ry));
  end

  function [pf, m] = idtlz1pfront(m)
    [pf, m] = utils.usimplex(nobj, m);
    pf = (1 - pf) / 2; 
  end
end



