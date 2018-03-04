function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = wfg9(sizes, k, l, doval)
% WFG8 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% Suggested nobj = 3, n = nobj + 9 = 12.
% 0 <= xi <= 2i.
% 1 <= k <= n. k is the position parameter and must be divisible by nobj - 1.
% 1 <= l <= n - k. l is the distance parameter. It does not need to be
% divisible by 2.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

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
t1 = @(y) [cell2mat(cellfun(@(i)... 
            mop.wfg.b_param(y, i, i + 1 : n, ones(1, n - i), 0.98 / 49.98, 0.02, 50),... 
            num2cell(1 : n - 1), 'UniformOutput', false)),...
           y(:, n)];

t2 = @(y) [mop.wfg.s_decept(y(:, 1 : k), 0.35, 0.001, 0.05),...
           mop.wfg.s_multi(y(:, k + 1 : n), 30, 95, 0.35)];

t3 = @(y) [cell2mat(cellfun(@(i)... 
            mop.wfg.r_nonsep(y, (i - 1) * k / (M - 1) + 1 : i * k / (M - 1), k / (M - 1)),... 
            num2cell(1 : M - 1), 'UniformOutput', false)),...
           mop.wfg.r_nonsep(y, k + 1 : n, l)];

t = {t1, t2, t3};

% wfg
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, ~, pfront] = mop.wfg.wfg(sizes, doval, h, t);
pset = @wfg9pset;
opts.MOPName = 'WFG9';

  function [ps, m] = wfg9pset(m)
    [ps, m] = utils.ugrid(M - 1, m);
    ps = [ps, zeros(m, n - M + 1)];
    
    ps(:, n) = 0.35;
    
    for i = n - 1 : - 1 : M
      u = mop.wfg.r_sum(ps, i + 1 : n, ones(1, n - i));
      v = 0.35.^((0.02 + 1.96 * u).^-1);
      ps(:, i) = v;
    end
    
    ps = bsxfun(@times, 2 * (1 : n), ps);
  end
end


