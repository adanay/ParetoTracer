% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = c3_dtlz4(sizes, alpha, doval)
% Constrained DTLZ4 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% alpha > 0.
% Suggested nobj = 3, n = nobj + 9 = 12.

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
nc = nobj;

sizes = struct(...
  'variables', n,...
  'objectives', nobj,...
  'ineqnonlin', nc);

if nargin < 3
  doval = false;
end

if nargin < 2
  alpha = 1;
end

if alpha <= 0
  alpha = 1;
end

[objfun, lb, ub, lincon, ~, multfun, sizes, opts] = mop.dtlz.dtlz(sizes, @dtlz4f);
pset = @dtlz4pset;
pfront =  @dtlz4pfront;
opts.MOPName = 'C3DTLZ4';

% nonlinear constraints
nonlcon = @dtlz4c;
  function [fx] = dtlz4f(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    y = x(:, 1 : nobj - 1).^alpha;
    z = x(:, nobj : n);
    
    gz = sum((z - 0.5).^2, 2);  

    hy = cumprod([ones(m, 1), cos(pi * y / 2)], 2);
    hy = fliplr(hy);

    ry = [ones(m, 1), sin(pi * y(:, end : -1 : 1) / 2)];
    
    fx = bsxfun(@times, 1 + gz, hy .* ry);
  end

  function [cx] = dtlz4c(x)
    x = val.valmopargin(x, [], [], n, nc, doval);
    
    fx = dtlz4f(x);
    cx = 1 - bsxfun(@plus, fx.^2 / 4, sum(fx.^2, 2)) + fx.^2;
  end

  function [ps, m] = dtlz4pset(m)
    q = n - nobj + 1;
    [ps, m] = utils.ugrid(nobj - 1, m);
    ps = [ps, zeros(m, q)];
    pc = dtlz4c(ps);
    i = bsxfun(@gt, pc, 0);
    i = any(i, 2);
    ps(i, :) = [];
  end
  
  function [pf, m] = dtlz4pfront(m)
    [pf, m] = utils.usimplex(nobj, m);
    pf2 = pf.^2;
    pf = bsxfun(@rdivide, pf, sqrt(sum(pf2, 2) - 3 / 4 * max(pf2, [], 2))); 
  end
end



