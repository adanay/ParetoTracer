function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = zdt(sizes, doval, f1, g, h, c)
% ZDT multiobjective optimization problem.
%
% n > 1.
% lb = 0.
% ub = 1.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

n = sizes.variables;
nobj = sizes.objectives;

% objective function
objfun = @zdtf;

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
opts = pt.defopts(n, 2);
opts.UseVectorized = true;

if nargin < 6
  c = 0;
end

% Pareto set/front
pset = @zdtpset;
pfront = @zdtpfront;

  function [fx] = zdtf(x)
    x = val.valmopargin(x, [], [], n, nobj, doval);
    
    f1x = f1(x);
    gx = g(x);
    f2x = gx .* h(f1x, gx);
    
    fx = [f1x, f2x];
  end

  function [ps, m] = zdtpset(m)
    ps = zeros(m, n); % each row is an individual
    ps(:, 1) = linspace(0, 1, m); % for the ps set all variables to zero except the first one. 
  end

  function [pf, m] = zdtpfront(m)
    pf = linspace(0, 1, m)';
    pf = [pf, 1 - pf.^c];
  end
end

