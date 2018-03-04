function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = dtlz8(sizes, doval)
% DTLZ8 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% Suggested nobj = 3, n = 10 * nobj = 30.

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
nc = nobj;

sizes = struct(...
  'variables', n,...
  'objectives', nobj,...
  'ineqnonlin', nc);

if nargin < 2
  doval = false;
end

[objfun, lb, ub, lincon, ~, multfun, sizes, opts] = mop.dtlz.dtlz(sizes, @dtlz8f);
pset = @dtlz8pset;
pfront = @dtlz8pfront;
opts.MOPName = 'DTLZ8';

% nonlinear constraints
nonlcon = @dtlz8c;

  function [fx] = dtlz8f(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    fx = zeros(m, nobj);
    for i = 1 : nobj
      fx(:, i) = mean(x(:, (i - 1) * n / nobj + 1 : i * n / nobj), 2);
    end
  end

  function [cx] = dtlz8c(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nc, doval);
    fx = dtlz8f(x);
    
    cx = zeros(m, nc);
    cx(:, 1 : nc - 1) = 1 - bsxfun(@plus, fx(:, nc), 4 * fx(:, 1 : nc - 1));
    
    if nc == 2
      cx(:, nc) = 0;
    else
      [~, rank] = sort(fx(:, 1 : nc - 1), 2);
      cx(:, nc) = 1 - 2 * fx(:, nc) - fx(rank(:, 1)) - fx(rank(:, 2));
    end
  end

  function [ps, m] = dtlz8pset(m)
    if nobj == 2
      ps = linspace(0, 1, m)';
      ps = [(1 - ps) / 4, ps];
     else
      temp = utils.usimplex(3, m / (nobj - 1));
      temp(:, 3) = temp(:, 3) / 2;
      temp = temp(temp(:, 1) >= (1 - temp(:, 3)) / 4 & temp(:, 1) <= temp(:, 2) & temp(:, 3) <= 1/3, :);

      ps = [repmat(temp(:, 2), nobj - 1, nobj - 1), repmat(temp(:, 3), nobj - 1, 1)];
      for i = 1 : nobj - 1
        ps((i - 1) * size(temp, 1) + 1 : i * size(temp, 1), i) = temp(:, 1);
      end

      gap = sort(unique(ps(:, nobj)));
      gap = gap(2) - gap(1);

      temp = (1 / 3 : gap : 1)';
      ps = [ps; repmat((1 - temp) / 4, 1, nobj - 1), temp];
      ps = unique(ps, 'rows');
    end
    m = size(ps, 1);
    ps = [ps, zeros(m, n - nobj)];
  end

  function [pf, m] = dtlz8pfront(m)
    if nobj == 2
      pf = linspace(0, 1, m)';
      pf = [(1 - pf) / 4, pf];
     else
      temp = utils.usimplex(3, m / (nobj - 1));
      temp(:, 3) = temp(:, 3) / 2;
      temp = temp(temp(:, 1) >= (1 - temp(:, 3)) / 4 & temp(:, 1) <= temp(:, 2) & temp(:, 3) <= 1/3, :);

      pf = [repmat(temp(:, 2), nobj - 1, nobj - 1), repmat(temp(:, 3), nobj - 1, 1)];
      for i = 1 : nobj - 1
        pf((i - 1) * size(temp, 1) + 1 : i * size(temp, 1), i) = temp(:, 1);
      end

      gap = sort(unique(pf(:, nobj)));
      gap = gap(2) - gap(1);

      temp = (1 / 3 : gap : 1)';
      pf = [pf; repmat((1 - temp) / 4, 1, nobj - 1), temp];
      pf = unique(pf, 'rows');
    end
    m = size(pf, 1);
  end
end



