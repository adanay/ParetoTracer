% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = wfg2(sizes, k, l, doval)
% WFG2 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% Suggested nobj = 3, n = nobj + 9 = 12.
% 0 <= xi <= 2i.
% 1 <= k <= n. k is the position parameter and must be divisible by nobj - 1.
% 1 <= l <= n - k. l is the distance parameter and should be divisible by 2.

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
  h = @(x) [mop.wfg.convex1(x(:, 1 : M - 1)), mop.wfg.discM(x(:, 1 : M - 1), 1, 5, 1)]; % (M - 1) => 2
else
  h = @(x) [mop.wfg.convex1(x(:, 1 : M - 1)), mop.wfg.convexI(x(:, 1 : M - 1)), mop.wfg.discM(x(:, 1 : M - 1), 1, 5, 1)]; % (M - 1) => M
end

% transformation functions
t1 = @(y) [y(:, 1 : k), mop.wfg.s_linear(y(:, k + 1 : n), 0.35)]; % n => n
t2 = @(y) [y(:, 1 : k),...
           cell2mat(cellfun(@(i)...
             mop.wfg.r_nonsep(y, k + 2 * (i - k) - 1 : k + 2 * (i - k), 2),... 
             num2cell(k + 1 : k + l / 2), 'UniformOutput', false))]; % n => ?

t3 = @(y) [cell2mat(cellfun(@(i)... 
            mop.wfg.r_sum(y, (i - 1) * k / (M - 1) + 1 : i * k / (M - 1), ones(1, k / (M - 1))),... 
            num2cell(1 : M - 1), 'UniformOutput', false)),...
           mop.wfg.r_sum(y, k + 1 : k + l / 2, ones(1, l / 2))]; % ? => M 

t = {t1, t2, t3};

% wfg
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts] = mop.wfg.wfg(sizes, doval, h, t);
pset = @wfg2pset;
pfront = @wfg2pfront;
opts.MOPName = 'WFG2';

  function [ps, m] = wfg2pset(m)
    [ps, m] = utils.ugrid(M - 1, m);
    ps(:, 1) = sqrt(ps(:, 1));
    ps = bsxfun(@times, 2 * (1 : (M - 1)), ps);
    ps(:, M : n) = repmat(2 * (M : n) * 0.35, m, 1);
    pf = objfun(ps);
    i = utils.ndsort(pf, 1) == 1;
    ps(~i, :) = []; 
  end

  function [pf, m] = wfg2pfront(m)
    [pf, m] = utils.ugrid(M - 1, m);
    pf(:, 1) = sqrt(pf(:, 1));
    temp = [pf, zeros(m, 1)];
    pf = mop.wfg.convexI(temp);
    pf = [pf, mop.wfg.discM(temp, 1, 5, 1)];
    pf = bsxfun(@times, 2 : 2 : 2 * M, pf);
    pf = pf(utils.ndsort(pf, 1) == 1, :);
  end
end


