% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = c2_dtlz2(sizes, doval)
% Constrained DTLZ2 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
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
nc = 1;

sizes = struct(...
  'variables', n,...
  'objectives', nobj,...
  'ineqnonlin', nc);

if nargin < 2
  doval = false;
end

[objfun, lb, ub, lincon, ~, multfun, sizes, opts] = mop.dtlz.dtlz(sizes, @dtlz2f);
pset = @dtlz2pset;
pfront = @dtlz2pfront;
opts.MOPName = 'C2DTLZ2';

% nonlinear constraints
nonlcon = @dtlz2c;

  function [fx] = dtlz2f(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    y = x(:, 1 : nobj - 1);
    z = x(:, nobj : n);
    
    gz = sum((z - 0.5).^2, 2);

    hy = cumprod([ones(m, 1), cos(pi * y / 2)], 2);
    hy = fliplr(hy);
    
    ry = [ones(m, 1), sin(pi * y(:, end : -1 : 1) / 2)];
    
    fx = bsxfun(@times, 1 + gz, hy .* ry);
  end

  function [cx] = dtlz2c(x)
    x = val.valmopargin(x, [], [], n, nc, doval);
    
    if nobj == 3 
      r = 0.4; 
    else
      r = 0.5;
    end
    
    fx = dtlz2f(x);
    cx = min(min(bsxfun(@plus, (fx - 1).^2, sum(fx.^2, 2)) - fx.^2 - r^2, [], 2), sum((fx - 1 / sqrt(nobj)).^2, 2) - r^2);
  end

  function [ps, m] = dtlz2pset(m)
    q = n - nobj + 1;
    [ps, m] = utils.ugrid(nobj - 1, m);
    ps = [ps, 0.5 * ones(m, q)]; 
    
    if nobj == 3 
      r = 0.4; 
    else
      r = 0.5;
    end
    
    pf = dtlz2f(ps);
    ps(min(min(bsxfun(@plus, (pf - 1).^2, sum(pf.^2, 2)) - pf.^2 - r^2, [], 2), sum((pf - 1 / sqrt(nobj)).^2, 2) - r^2) > 0, :) = [];
    m = size(ps, 1);
  end

  function [pf, m] = dtlz2pfront(m)
    pf = utils.usimplex(nobj, m);
    pf = bsxfun(@rdivide, pf, sqrt(sum(pf.^2, 2))); 
     
    if nobj == 3 
      r = 0.4; 
    else
      r = 0.5;
    end
     
    pf(min(min(bsxfun(@plus, (pf - 1).^2, sum(pf.^2, 2)) - pf.^2 - r^2, [], 2), sum((pf - 1 / sqrt(nobj)).^2, 2) - r^2) > 0, :) = [];
    m = size(pf, 1);
  end
end



