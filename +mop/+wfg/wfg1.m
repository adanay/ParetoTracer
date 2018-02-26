% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = wfg1(sizes, k, doval)
% WFG1 multiobjective optimization problem.
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
  h = @(x) [mop.wfg.convex1(x(:, 1 : M - 1)), mop.wfg.mixedM(x(:, 1 : M - 1), 1, 5)]; % (M - 1) => 2
else
  h = @(x) [mop.wfg.convex1(x(:, 1 : M - 1)), mop.wfg.convexI(x(:, 1 : M - 1)), mop.wfg.mixedM(x(:, 1 : M - 1), 1, 5)]; % (M - 1) => M
end

% transformation functions
t1 = @(y) [y(:, 1 : k), mop.wfg.s_linear(y(:, k + 1 : n), 0.35)]; % n => n
t2 = @(y) [y(:, 1 : k), mop.wfg.b_flat(y(:, k + 1 : n), 0.8, 0.75, 0.85)]; % n => n
t3 = @(y) mop.wfg.b_poly(y, 0.02); % n => n
t4 = @(y) [cell2mat(cellfun(@(i)... 
             mop.wfg.r_sum(y, (i - 1) * k / (M - 1) + 1 : i * k / (M - 1), 2 * ((i - 1) * k / (M - 1) + 1 : i * k / (M - 1))),... 
             num2cell(1 : M - 1), 'UniformOutput', false)),...
           mop.wfg.r_sum(y, k + 1 : n, 2 * (k + 1 : n))]; % n => M

t = {t1, t2, t3, t4};

% wfg
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts] = mop.wfg.wfg(sizes, doval, h, t);
pset = @wfg1pset;
pfront = @wfg1pfront;
opts.MOPName = 'WFG1';

  function [ps, m] = wfg1pset(m)
    [ps, m] = utils.ugrid(M - 1, m);
    ps = ps.^(50);
    ps(:, M : n) = repmat(2 * (M : n) * 0.35, m, 1);
  end

  function [pf, m] = wfg1pfront(m)
    [pf, m] = utils.usimplex(M, m);
    for i = 1 : m
      c = ones(1, M);
      p = find(pf(i, :) ~= 0,1);
      for j = p + 1 : M
        temp = pf(i, j) / pf(i, p) * prod(1 - c(M - j + 2 : M - p));
        c(M - j + 1) = (temp^2 - temp + sqrt(2 * temp)) / (temp^2 + 1);
      end
      for j = 1 : M
        pf(i, j) = prod(1 - c(1 : M - j)) .* (1 - sqrt(1 - c(M - j + 1)^2));
      end
      temp = acos(c(1)) * 2 / pi;                   
      pf(i, M) = 1 - temp - cos(10 * pi * temp + pi / 2) / 10 / pi;
    end
    pf = bsxfun(@times, 2 : 2 : 2 * M, pf);
  end
end


