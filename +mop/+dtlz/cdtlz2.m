% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = cdtlz2(sizes, doval)
% Convex DTLZ2 multiobjective optimization problem.
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

if nargin < 2
  doval = false;
end

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset] = mop.dtlz.dtlz(sizes, @cdtlz2f);
pfront = @cdtlz2pfront;
opts.MOPName = 'CDTLZ2';

  function [fx] = cdtlz2f(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    y = x(:, 1 : nobj - 1);
    z = x(:, nobj : n);
    
    gz = sum((z - 0.5).^2, 2);

    hy = cumprod([ones(m, 1), cos(pi * y / 2)], 2);
    hy = fliplr(hy);
    
    ry = [ones(m, 1), sin(pi * y(:, end : -1 : 1) / 2)];
    
    fx = bsxfun(@times, 1 + gz, hy .* ry);
    fx = [fx(:, 1 : nobj - 1).^4, fx(:, nobj).^2];
  end

  function [pf, m] = cdtlz2pfront(m)
    [pf, m] = utils.usimplex(nobj, m);
    pf = pf.^2;
    temp = sum(sqrt(pf(:, 1 : end - 1)), 2) + pf(:, end);
    pf = pf ./ [repmat(temp.^2, 1, size(pf, 2) - 1), temp];
  end
end



