function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = zdt4(sizes, doval)
% ZDT4 multiobjective optimization problem.
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

f1 = @(x) x(:, 1);
g = @(x) 1 + 10 * (n - 1) + sum(x(:, 2 : n).^2 - 10 * cos(4 * pi * x(:, 2 : n)), 2);
h = @(f1x, gx) 1 - sqrt(f1x ./ gx);

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.zdt.zdt(sizes, doval, f1, g, h, 0.5);
opts.MOPName = 'ZDT4';
end