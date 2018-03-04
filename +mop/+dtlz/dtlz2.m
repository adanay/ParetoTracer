function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = dtlz2(sizes, doval)
% DTLZ2 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% Suggested nobj = 3, n = nobj + 9 = 12.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

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

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset] = mop.dtlz.dtlz(sizes, @dtlz2f);
pfront = @dtlz2pfront;
opts.MOPName = 'DTLZ2';

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

  function [pf, m] = dtlz2pfront(m)
    [pf, m] = utils.usimplex(nobj, m);
    pf = bsxfun(@rdivide, pf, sqrt(sum(pf.^2, 2))); 
  end
end



