% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = wfg6(sizes, k, l, doval)
% WFG6 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% Suggested nobj = 3, n = nobj + 9 = 12.
% 0 <= xi <= 2i.
% 1 <= k <= n. k is the position parameter and must be divisible by nobj - 1.
% 1 <= l <= n - k. l is the distance parameter. It does not need to be
% divisible by 2.

% problem sizes
[n, nobj] = val.valmopsizes(sizes, true);
if n < 3
  n = 3;
end
if nobj < 2
  nobj = 2;
end
if nobj > n - 1
  nobj = n - 1;
end

sizes = struct(...
  'variables', n,...
  'objectives', nobj);

if nargin < 4
  doval = false;
end

if nargin < 2
  k = floor(n / 2);
  k = k - mod(k, nobj - 1);
  if k < nobj - 1
    k = nobj - 1;
  end
end

if nargin < 3
  l = n - k;
  l = l - mod(l, 2);
end

M = nobj;

% shape functions
if M == 2
  h = @(x) [mop.wfg.concave1(x(:, 1 : M - 1)), mop.wfg.concaveM(x(:, 1 : M - 1))]; % (M - 1) => 2
else
  h = @(x) [mop.wfg.concave1(x(:, 1 : M - 1)), mop.wfg.concaveI(x(:, 1 : M - 1)), mop.wfg.concaveM(x(:, 1 : M - 1))]; % (M - 1) => M
end

% transformation functions
t1 = @(y) [y(:, 1 : k), mop.wfg.s_linear(y(:, k + 1 : n), 0.35)];
t2 = @(y) [cell2mat(cellfun(@(i)... 
             mop.wfg.r_nonsep(y, (i - 1) * k / (M - 1) + 1 : i * k / (M - 1), k / (M - 1)),... 
             num2cell(1 : M - 1), 'UniformOutput', false)),...
           mop.wfg.r_nonsep(y, k + 1 : n, l)];

t = {t1, t2};

% wfg
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.wfg.wfg(sizes, doval, h, t);
opts.MOPName = 'WFG6';
end


