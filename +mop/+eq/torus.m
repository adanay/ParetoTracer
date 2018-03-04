function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = torus(sizes, doval, a, r, R)
% Quadratic problem with one torus equality constraint.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if nargin < 2
  doval = false;
end

if nargin < 3
  [n, nobj] = val.valmopsizes(sizes, true);
  if (n < 3)
    n = 3;
  end
  if n < nobj
    n = nobj;
  end
  
  if nobj == 2
    a = [ones(n, 1), -ones(n, 1)];
  elseif nobj == 3
    a = ones(n, 1);
    a(2 : 2 : end) = -1;
    a = [ones(n, 1), -ones(n, 1), a];
  else
    error('mop:eq:quad1', 'The vector a of size (n x nobj) should be specified.');
  end
else
  a = a(:, :);  
  n = size(a, 1);
  nobj = size(a, 2);
  
  if n < nobj
    a = a';
    n = size(a, 1);
    nobj = size(a, 2);
  end
end
a = shiftdim(a, -1);

if nargin < 4
  r = 0.3;
end

if nargin < 5
  R = 0.5;
end

sizes = struct(...
  'variables', n,... 
  'objectives', nobj,...
  'ineqlin', 0,...
  'eqlin', 0,...
  'ineqnonlin', 0,...
  'eqnonlin', 1);

% objective function, Jacobian, and Hessians
objfun = struct('f', @torusf, 'J', @torusJ, 'H', @torusH);

% box constraints
lb = [];
ub = [];

% linear constraints
lincon = [];

% nonlinear constraints
nonlcon = struct('ceq', @torusceq, 'Jceq', @torusJceq, 'Hceq', @torusHceq);                

% multiply functions
multfun = [];

% options
opts = pt.defopts(n, nobj);
opts.UseVectorized = true;
opts.MOPName = 'Torus';

% Pareto set and front
pset = [];
pfront = [];

  function [fx] = torusf(x)
    x = val.valmopargin(x, [], [], n, nobj, doval);
    
    fx = sum(bsxfun(@minus, x, a).^2, 2);
    fx = permute(fx, [1 3 2]);
  end

  function [Jx] = torusJ(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    Jx = 2 * bsxfun(@minus, x, a);
    Jx = permute(Jx, [1 3 2]);
    
    if m == 1
      Jx = shiftdim(Jx, 1); % (nobj x n)
    end
  end

  function [Hx] = torusH(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    Hx = 2 * utils.eyen([m, n, n, nobj], 2, 3);
         
    if m == 1
      Hx = shiftdim(Hx, 1); % (n x n x nobj)
    end
  end

  function [ceqx] = torusceq(x)
    x = val.valmopargin(x, [], [], n, 1, doval);
    
    ceqx = r^2 - x(:, 3).^2 - (sqrt(x(:, 1).^2 + x(:, 2).^2) - R).^2;
  end

  function [Jceqx] = torusJceq(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, 1, doval);
    
    q = sqrt(x(:, 1).^2 + x(:, 2).^2);
    
    Jceqx = [-2 * x(:, 1) .* (q - R) ./ q,...
             -2 * x(:, 2) .* (q - R) ./ q,...
             -2 * x(:, 3),...
             zeros(m, n - 3)];
    
    if m > 1
      Jceqx = permute(Jceqx, [1 3 2]);
    end
  end

  function [Hceqx] = torusHceq(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, 1, doval);
    
    q = sqrt(x(:, 1).^2 + x(:, 2).^2);
    
    r1 = [-2 + 2 * R ./ q - 2 * x(:, 1).^2 * R ./ q.^3, -2 * x(:, 1) .* x(:, 2) * R ./ q.^3, zeros(m, n - 2)];
    r2 = [-2 * x(:, 1) .* x(:, 2) * R ./ q.^3, -2 + 2 * R ./ q - 2 * x(:, 2).^2 * R ./ q.^3, zeros(m, n - 2)];
    r3 = [zeros(m, 2), -2 * ones(m, 1), zeros(m, n - 3)];
    
    Hceqx = cat(2, permute(r1, [1 3 2]),...
                   permute(r2, [1 3 2]),...
                   permute(r3, [1 3 2]),...
                   zeros(m, n - 3, n));
    
    if m == 1
      Hceqx = shiftdim(Hceqx, 1); % (n x n x nceq)
    end
  end
end

