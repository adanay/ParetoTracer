% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset] = dtlz(sizes, dtlzf, c)
% DTLZ multiobjective optimization problem.
%
% n > 1.
% lb = 0.
% ub = 1.

n = sizes.variables;
nobj = sizes.objectives;

% objective function
objfun = dtlzf;

% box constraints
lb = zeros(n, 1);
ub = ones(n, 1);

% linear constraints
lincon = [];

% nonlinear constraints
nonlcon = [];

% multiply functions
multfun = [];

% def optimization opts
opts = pt.defopts(n, nobj);
opts.UseVectorized = true;

% Pareto set/front
pset = @dtlzpset;

if nargin < 3
  c = 0.5;
end

  function [ps, m] = dtlzpset(m)
    q = n - nobj + 1;
    [ps, m] = utils.ugrid(nobj - 1, m);
    ps = [ps, c * ones(m, q)];
  end
end

