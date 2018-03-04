function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = idtlz2(sizes, doval)
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

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset] = mop.dtlz.dtlz(sizes, @idtlz2f);
pfront = @idtlz2pfront;
opts.MOPName = 'IDTLZ2';

  function [fx] = idtlz2f(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    y = x(:, 1 : nobj - 1);
    z = x(:, nobj : n);
    
    gz = sum((z - 0.5).^2, 2);

    hy = cumprod([ones(m, 1), cos(pi * y / 2)], 2);
    hy = fliplr(hy);
    
    ry = [ones(m, 1), sin(pi * y(:, end : -1 : 1) / 2)];
    
    fx = bsxfun(@minus, 1 + gz, bsxfun(@times, 1 + gz, hy .* ry));
  end

  function [pf, m] = idtlz2pfront(m)
    [pf, m] = utils.usimplex(nobj, m);
    pf = bsxfun(@rdivide, pf, sqrt(sum(pf.^2, 2))); 
    pf = 1 - pf;
  end
end



