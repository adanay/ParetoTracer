function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = dd2(sizes, doval)
% NBI test problem modification.
% n = 5, nobj = 2
% two equalities, -2 <= x <= 2

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if nargin < 2
  doval = false;
end

n = 5;
nobj = 2;

sizes = struct(...
  'variables', n,... 
  'objectives', nobj,...
  'ineqlin', 0,...
  'eqlin', 1,...
  'ineqnonlin', 0,...
  'eqnonlin', 1);
 
% objective function, Jacobian, and Hessians
objfun = struct('f', @dd2f, 'J', @dd2J, 'H', @dd2H);

% box constraints
lb = -2 * ones(n, 1);
ub = 2 * ones(n, 1);

% linear constraints
lincon = struct(...
'A', [],...
'b', [],...
'Aeq', [1, 2, -1, -0.5, 1],...
'beq', -2);

% nonlinear constraints
nonlcon = struct('ceq', @dd2ceq, 'Jceq', @dd2Jceq, 'Hceq', @dd2Hceq);  
               
multfun = [];

% options
opts = pt.defopts(n, nobj);
opts.UseVectorized = true;
opts.MOPName = 'DD2';

% Pareto set and front
pset = [];
pfront = [];

  function [fx] = dd2f(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    f1x = x(:, 1).^2 + x(:, 2).^2 + x(:, 3).^2 + x(:, 4).^2 + x(:, 5).^2;
    f2x = 3 * x(:, 1) + 2 * x(:, 2) - x(:, 3) / 3 + 0.01 * (x(:, 4) - x(:, 5)).^3;
    
    if m > 1
      fx = [f1x, f2x]; % (m x 2)
    else
      fx = [f1x; f2x]; % (2 x 1)
    end
  end

  function [Jx] = dd2J(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    df1x = 2 * x;
    d = repmat([3, 2, -1 / 3], m, 1);
    df2x = [d, 0.03 * (x(:, 4) - x(:, 5)).^2, -0.03 * (x(:, 4) - x(:, 5)).^2]; 
    
    Jx = cat(3, df1x, df2x);
    Jx = permute(Jx, [1, 3, 2]); % (m x nobj x n)
    if m == 1
      Jx = shiftdim(Jx, 1); % (nobj x n)
    end
  end

  function [Hx] = dd2H(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    d =  0.06 * (x(:, 4) - x(:, 5));
    
    ddf1x = repmat(shiftdim(2 * eye(n), -1), m, 1); % (m x n x n)
    ddf2x = cat(2, zeros(m, 3, n),... 
                   cat(3, zeros(m, 1, 3), d, -d),...
                   cat(3, zeros(m, 1, 3), -d, d));
         
    Hx = cat(4, ddf1x, ddf2x); % (m x n x n x nobj)
    if m == 1
      Hx = shiftdim(Hx, 1); % (n x n x nobj)
    end
  end

  function [ceqx] = dd2ceq(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, 1, doval);
    
    ceqx = 4 * x(:, 1) - 2 * x(:, 2) + 0.8 * x(:, 3) + 0.6 * x(:, 4) + 0.5 * x(:, 5).^2;
    
    if m == 1
      ceqx = ceqx'; % (nceq x 1)
    end
  end

  function [Jceqx] = dd2Jceq(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, 1, doval);
    
    Jceqx = [repmat([4, -2, 0.8, 0.6], m, 1), x(:, 5)]; % (m x 1 x n)
    
    if m > 1
      Jceqx = permute(Jceqx, [1 3 2]); % (m x 1 x n)
    end
  end

  function [Hceqx] = dd2Hceq(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, 1, doval);
    
    Hceqx = [zeros(n - 1, n); [zeros(1, n - 1), 1]]; % (n x n)
    
    if m ~= 1
      Hceqx = repmat(shiftdim(Hceqx, -1), m, 1); % (m x n x n)
    end
  end
end