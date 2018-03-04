function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = zdt6(sizes, doval)
% ZDT6 multiobjective optimization problem.
% Suggested nobj = 2, n = 10.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% problem sizes
n = val.valmopsizes(sizes, true);
if n < 2
  n = 2;
end
nobj = 2;
sizes = struct(...
'variables', n,... 
'objectives', nobj);

if nargin < 2
  doval = false;
end

f1 = @(x) 1 - exp(-4 * x(:, 1)) .* sin(6 * pi * x(:, 1)).^6;
g = @(x) 1 + 9 * (sum(x(:, 2 : n), 2) / (n - 1)).^0.25;
h = @(f1x, gx) 1 - (f1x ./ gx).^2;

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset] = mop.zdt.zdt(sizes, doval, f1, g, h);
pfront = @zdt6pfront;
opts.MOPName = 'ZDT6';

  function [pf, m] = zdt6pfront(m)
    pf = linspace(0.280775, 1, m)';
    pf = [pf, 1 - pf.^2];
  end
end