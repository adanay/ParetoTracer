% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = dtlz9(sizes, doval)
% DTLZ9 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% Suggested nobj = 3, n = 10 * nobj = 30.

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
nc = nobj - 1;

sizes = struct(...
  'variables', n,...
  'objectives', nobj,...
  'ineqnonlin', nc);

if nargin < 2
  doval = false;
end

[objfun, lb, ub, lincon, ~, multfun, sizes, opts] = mop.dtlz.dtlz(sizes, @dtlz9f);
pset = @dtlz9pset;
pfront = @dtlz9pfront;
opts.MOPName = 'DTLZ9';

% nonlinear constraints
nonlcon = @dtlz9c;

  function [fx] = dtlz9f(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    x = x.^0.1;
    
    fx = zeros(m, nobj);
    for i = 1 : nobj
      fx(:, i) = sum(x(:, floor((i - 1) * n / nobj) + 1 : floor(i * n / nobj)), 2);
    end
  end

  function [cx] = dtlz9c(x)
    x = val.valmopargin(x, [], [], n, nc, doval);
    fx = dtlz9f(x);
    
    cx = 1 - bsxfun(@plus, fx(:, nobj).^2, fx(:, 1 : nobj - 1).^2);
  end

  function [ps, m] = dtlz9pset(m)
    ps = linspace(0, 1, m)';
    ps = [repmat(cos(0.5 * pi * ps), 1, nobj - 1), sin(0.5 * pi * ps)];
    ps = ps.^10;
    ps = [ps, zeros(m, n - nobj)];
  end

  function [pf, m] = dtlz9pfront(m)
    pf = linspace(0, 1, m)';
    pf = [repmat(cos(0.5 * pi * pf), 1, nobj - 1), sin(0.5 * pi * pf)];
  end
end



