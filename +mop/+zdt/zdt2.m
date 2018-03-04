function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = zdt2(sizes, doval)
% ZDT2 multiobjective optimization problem.
% Suggested nobj = 2, n = 30.

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

f1 = @(x) x(:, 1);
g = @(x) 1 + 9 * sum(x(:, 2 : n), 2) / (n - 1);
h = @(f1x, gx) 1 - (f1x ./ gx).^2;

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.zdt.zdt(sizes, doval, f1, g, h, 2);
opts.MOPName = 'ZDT2';
end
