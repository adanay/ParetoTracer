% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = wfg(sizes, doval, h, t, A, S, D, zmax)
% General WFG function.
%
% n is the number of variables.
% M <= n. M is the number of objectives.
% h represents the shape functions. It is a handle to a function that
% receives M - 1 parameters and returns M values.
% t is a cell array of p handles to the transformation functions. The
% functions will be applied this way: tp(...(t3(t2(t1(z))). The output of 
% any function must be compatible with the input of the function to its 
% right. t1 must accept n parameters and tp must return M values.
%
% A in {0, 1} is a vector of M - 1 degeneracy constants. For each A(i) = 0
% the dimensionality of the Pareto front is reduced by one. The default
% value is the vector of all ones.
% S is a vector of M scaling constants. The default value is 2 * (1 : M)'.  
% D > 0 is a distant scaling constant. The default value is one.
% zmax is the vector of the upper bounds for the variables. The lower
% bounds are assumed to be zero. The default value is 2 * (1 : n)'.

n = sizes.variables;
M = sizes.objectives; % nobj

% default values
if ~exist('A', 'var')
  A = ones(1, M - 1);
else
  A = A(:)';
end
if ~exist('S', 'var')
  S = 2 * (1 : M);
else
  S = S(:)';
end
if ~exist('D', 'var')
  D = 1;
end
if ~exist('zmax', 'var')
  zmax = 2 * (1 : n);
end

z01 = @(z) bsxfun(@rdivide, z, zmax); % normalizing function
tp = @(z01) pipeline(z01, t); % transformation function

% underlying parameters
x = @(tp) [bsxfun(@max, tp(:, M), A) .* (tp(:, 1 : M - 1) - 0.5) + 0.5, tp(:, M)]; 

% objective function
f = @(x) bsxfun(@plus, bsxfun(@times, S, h(x(:, 1 : M - 1))), D * x(:, M)); 

% generic WFG function.
objfun = @wfgf;

% box constraints
lb = zeros(n, 1);
ub = zmax;

% linear constraints
lincon = [];

% nonlinear constraints
nonlcon = [];

% multiply functions
multfun = [];

% def optimization opts
opts = pt.defopts(n, M);
opts.UseVectorized = true;

% Pareto set/front
pset = @wfgpset;
pfront = @wfgpfront;

  function [fz] = wfgf(z)
    z = val.valmopargin(z, [], [], n, M, doval);
    fz = pipeline(z, {z01, tp, x, f}); % @(z) f(x(tp(z01(z))));
  end

  function [ps, m] = wfgpset(m)
    [ps, m] = utils.ugrid(M - 1, m);
    ps = bsxfun(@times, 2 * (1 : M - 1), ps);
    ps(:, M : n) = repmat(2 * (M : n) * 0.35, m, 1);
  end

  function [pf, m] = wfgpfront(m)
    [pf, m] = utils.usimplex(M, m);
    pf = bsxfun(@rdivide, pf, sqrt(sum(pf.^2, 2)));
    pf = bsxfun(@times, 2 : 2 : 2 * M, pf);
  end
end


function [tp] = pipeline(z01, t)
tp = z01;
p = length(t);
for j = 1 : p
  tp = t{j}(tp);
end
end