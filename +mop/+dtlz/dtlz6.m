% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = dtlz6(sizes, k, doval)
% DTLZ6 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% k is the position parameter. 2 <= k <= nobj.
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
  k = 1;
end

if k < 2
  k = 2;
end

if k > nobj
  k = nobj;
end

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset] = mop.dtlz.dtlz(sizes, @dtlz6f, 0);
pfront = @dtlz6pfront;
opts.MOPName = 'DTLZ6';

  function [fx] = dtlz6f(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    y = x(:, 1 : nobj - 1);
    z = x(:, nobj : n);
    
    gz = sum(z.^0.1, 2);
    
    a = [y(:, 1 : k - 1), bsxfun(@rdivide, (1 + 2 * bsxfun(@times, gz, y(:, k : end))), (2 * (1 + gz)))]; 
    
    hy = cumprod([ones(m, 1), cos(pi * a / 2)], 2);
    hy = fliplr(hy);
    
    ry = [ones(m, 1), sin(pi * a(:, end : -1 : 1) / 2)];
    
    fx = bsxfun(@times, 1 + gz, hy .* ry);
  end

  function [pf, m] = dtlz6pfront(m)
    pf = utils.ugrid(k, m);
    pf = bsxfun(@rdivide, pf, sqrt(sum(pf.^2, 2)));
    pf = [pf(:, ones(1, nobj - k)), pf];
    pf = bsxfun(@rdivide, pf, sqrt(2).^max([nobj - k, nobj - k : -1 : 2 - k], 0));
  end
end



