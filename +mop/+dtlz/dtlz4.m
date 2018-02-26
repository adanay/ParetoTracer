% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = dtlz4(sizes, alpha, doval)
% DTLZ4 multiobjective optimization problem.
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

sizes = struct(...
  'variables', n,...
  'objectives', nobj);

if nargin < 3
  doval = false;
end

if nargin < 2
  alpha = 1;
end

if alpha <= 0
  alpha = 1;
end

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset] = mop.dtlz.dtlz(sizes, @dtlz4f);
pfront = @dtlz4pfront;
opts.MOPName = 'DTLZ4';

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

  function [pf, m] = dtlz4pfront(m)
    [pf, m] = utils.usimplex(nobj, m);
    pf = bsxfun(@rdivide, pf, sqrt(sum(pf.^2, 2))); 
  end
end



