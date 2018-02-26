% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = quad2(sizes, doval, a, c, r)
% Quadratic problem with one quadratic equality constraint.

if nargin < 2
  doval = false;
end

if nargin < 3
  [n, nobj] = val.valmopsizes(sizes, true);
  if n < nobj
    n = nobj;
  end
  
  if nobj == 2
    n1 = floor(n / 2);
    n2 = n - n1;
    a = [ones(n, 1), [-ones(n1, 1); ones(n2, 1)]] / sqrt(n);
  else
    error('mop:eq:quad2', 'The vector a of size (n x nobj) should be specified.');
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
  c = zeros(n, 1);
end
c = c(:)';

if nargin < 5
  r = 1;
end

sizes = struct(...
  'variables', n,... 
  'objectives', nobj,...
  'ineqlin', 0,...
  'eqlin', 0,...
  'ineqnonlin', 0,...
  'eqnonlin', 1);

% objective function, Jacobian, and Hessians
objfun = struct('f', @quad2f, 'J', @quad2J, 'H', @quad2H);

% box constraints
lb = [];
ub = [];

% linear constraints
lincon = [];

% nonlinear constraints
nonlcon = struct('ceq', @quad2ceq, 'Jceq', @quad2Jceq, 'Hceq', @quad2Hceq);                
% multiply functions
multfun = [];

% options
opts = pt.defopts(n, nobj);
opts.UseVectorized = true;
opts.MOPName = 'Quad prob with nonlin eq';

% Pareto set and front
pset = [];
pfront = [];

  function [fx] = quad2f(x)
    x = val.valmopargin(x, [], [], n, nobj, doval);
    
    fx = sum(bsxfun(@minus, x, a).^2, 2);
    fx = permute(fx, [1 3 2]);
  end

  function [Jx] = quad2J(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    Jx = 2 * bsxfun(@minus, x, a);
    Jx = permute(Jx, [1 3 2]);
    
    if m == 1
      Jx = shiftdim(Jx, 1); % (nobj x n)
    end
  end

  function [Hx] = quad2H(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    Hx = 2 * utils.eyen([m, n, n, nobj], 2, 3);
         
    if m == 1
      Hx = shiftdim(Hx, 1); % (n x n x nobj)
    end
  end

  function [ceqx] = quad2ceq(x)
    x = val.valmopargin(x, [], [], n, 1, doval);
    
    ceqx = sum(bsxfun(@minus, x, c).^2, 2) - r^2;
    ceqx = permute(ceqx, [1 3 2]);
  end

  function [Jceqx] = quad2Jceq(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, 1, doval);
    
    Jceqx = 2 * bsxfun(@minus, x, c);
    Jceqx = permute(Jceqx, [1 3 2]);
    
    if m == 1
      Jceqx = shiftdim(Jceqx, 1); % (1 x n)
    end
  end

  function [Hceqx] = quad2Hceq(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, 1, doval);
    
    Hceqx = 2 * utils.eyen([m, n, n], 2, 3);
         
    if m == 1
      Hceqx = shiftdim(Hceqx, 1); % (n x n x nobj)
    end
  end
end

