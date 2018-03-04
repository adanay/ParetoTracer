function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = dent(sizes, doval, b)
% Quadratic problem.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if nargin < 2
  doval = false;
end
if nargin < 3
  b = 2;
end

lambda = 0.85;

n = 2;
nobj = 2;

sizes = struct(...
  'variables', n,...
  'objectives', nobj);

% objective function, Jacobian, and Hessians
objfun = struct('f', @dentf, 'J', @dentJ, 'H', @dentH);

% box constraints
lb = -[b, b];
ub = [b, b];

% linear constraints
lincon = [];

% nonlinear constraints
nonlcon = [];          

% multiply functions
multfun = [];

% options
opts = pt.defopts(n, nobj);
opts.UseVectorized = true;
opts.MOPName = 'Dent';

% Pareto set and front
pset = @dentpset;
pfront = [];

  function [fx] = dentf(x)
    x = val.valmopargin(x, [], [], n, nobj, doval);
    
    s1 = x(:, 1) + x(:, 2);
    s2 = x(:, 1) - x(:, 2); 

    q1 = sqrt(1 + s1.^2);
    q2 = sqrt(1 + s2.^2);

    e = 2 * lambda * exp(-s2.^2); 
    
    f1 = 0.5 * (q1 + q2 + s2 + e);
    f2 = 0.5 * (q1 + q2 - s2 + e);
    
    fx = [f1, f2];
  end

  function [Jx] = dentJ(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    s1 = x(:, 1) + x(:, 2);
    s2 = x(:, 1) - x(:, 2); 

    q1 = sqrt(1 + s1.^2);
    q2 = sqrt(1 + s2.^2);
    
    r1 = s1 ./ q1;
    r2 = s2 ./ q2;
    
    e = 2 * lambda * exp(-s2.^2);
    e1 = s2 .* e;
    
    df1 = [0.5 * (r1 + r2 + 1) - e1...
           0.5 * (r1 - r2 - 1) + e1];
    df2 = [0.5 * (r1 + r2 - 1) - e1...
           0.5 * (r1 - r2 + 1) + e1];    

    if m == 1
      Jx = [df1; df2];
    else
      Jx = cat(2, permute(df1, [1 3 2]), permute(df2, [1 3 2]));
    end
  end

  function [Hx] = dentH(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    s1 = x(:, 1) + x(:, 2);
    s2 = x(:, 1) - x(:, 2); 

    q1 = sqrt(1 + s1.^2);
    q2 = sqrt(1 + s2.^2);
    
    e = 2 * lambda * exp(-s2.^2);
    e2 = 2 * s2.^2 .* e;
    
    ddf1 = [0.5 * (1 ./ q1.^3 + 1 ./ q2.^3) - e + e2...
            0.5 * (1 ./ q1.^3 - 1 ./ q2.^3) + e - e2];
          
    ddf2 = [0.5 * (1 ./ q1.^3 - 1 ./ q2.^3) + e - e2...
            0.5 * (1 ./ q1.^3 + 1 ./ q2.^3) - e + e2];
    
    if m == 1
      Hx = [ddf1; ddf2];
      Hx = cat(3, Hx, Hx);
    else
      Hx = cat(2, permute(ddf1, [1 3 2]), permute(ddf2, [1 3 2]));
      Hx = cat(4, Hx, Hx);
    end
  end

  function [ps, m] = dentpset(m)
    ps = [linspace(lb(1), ub(1), m);...
          linspace(ub(2), lb(2), m)]';
  end
end

