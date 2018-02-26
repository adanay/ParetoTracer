% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = wfg7(sizes, k, doval)
% WFG7 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% Suggested nobj = 3, n = nobj + 9 = 12.
% 0 <= xi <= 2i.
% 1 <= k <= n. k is the position parameter and must be divisible by nobj - 1.

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
  k = floor(n / 2);
  k = k - mod(k, nobj - 1);
  if k < nobj - 1
    k = nobj - 1;
  end
end

M = nobj;

% shape functions
if M == 2
  h = @(x) [mop.wfg.concave1(x(:, 1 : M - 1)), mop.wfg.concaveM(x(:, 1 : M - 1))]; % (M - 1) => 2
else
  h = @(x) [mop.wfg.concave1(x(:, 1 : M - 1)), mop.wfg.concaveI(x(:, 1 : M - 1)), mop.wfg.concaveM(x(:, 1 : M - 1))]; % (M - 1) => M
end

% transformation functions
t1 = @(y) [cell2mat(cellfun(@(i)... 
             mop.wfg.b_param(y, i, i + 1 : n, ones(1, n - i), 0.98 / 49.98, 0.02, 50),... 
             num2cell(1 : k), 'UniformOutput', false)),...
           y(:, k + 1 : n)];
t2 = @(y) [y(:, 1 : k), mop.wfg.s_linear(y(:, k + 1 : n), 0.35)];
t3 = @(y) [cell2mat(cellfun(@(i)... 
             mop.wfg.r_sum(y, (i - 1) * k / (M - 1) + 1 : i * k / (M - 1), ones(1, k / (M - 1))),... 
             num2cell(1 : M - 1), 'UniformOutput', false)),...
           mop.wfg.r_sum(y, k + 1 : n, ones(1, n - k))];

t = {t1, t2, t3};

% wfg
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.wfg.wfg(sizes, doval, h, t);
opts.MOPName = 'WFG7';
end


