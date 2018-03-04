function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = sproblem(sizes, doval, a)
% S problem with one linear equality constraint.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if nargin < 2
  doval = false;
end

if nargin < 3
  [n, nobj] = val.valmopsizes(sizes, true);
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
    error('mop:misc:sproblem', 'The vector a of size (n x nobj) should be specified.');
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
diaga = diag(a)';
a = shiftdim(a, -1);

sizes = struct(...
  'variables', n,...
  'objectives', nobj);

% objective function, Jacobian, and Hessians
objfun = struct('f', @sproblemf, 'J', @sproblemJ, 'H', @sproblemH);

% box constraints
lb = [];
ub = [];

% linear constraints
lincon = [];

% nonlinear constraints
nonlcon = [];          

% multiply functions
multfun = [];

% options
opts = pt.defopts(n, nobj);
opts.UseVectorized = true;
opts.MOPName = 'SProblem';

% Pareto set and front
pset = [];
pfront = [];

  function [fx] = sproblemf(x)
    x = val.valmopargin(x, [], [], n, nobj, doval);
    
    trunx = x(:, 1 : nobj);
    trunfx = bsxfun(@minus, trunx, diaga).^4 - bsxfun(@minus, trunx, diaga).^2;
    
    fx = sum(bsxfun(@minus, x, a).^2, 2);
    fx = trunfx + permute(fx, [1 3 2]);
  end

  function [Jx] = sproblemJ(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    trunx = x(:, 1 : nobj);
    trunJx = 4 * bsxfun(@minus, trunx, diaga).^3 - 2 * bsxfun(@minus, trunx, diaga);
    
    Jx(m, n, nobj) = 0;
    for i = 1 : nobj
      Jx(:, i, i) = trunJx(:, i);
    end
    
    Jx = Jx + 2 * bsxfun(@minus, x, a);
    Jx = permute(Jx, [1 3 2]);
    
    if m == 1
      Jx = shiftdim(Jx, 1); % (nobj x n)
    end
  end

  function [Hx] = sproblemH(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    trunx = x(:, 1 : nobj);
    trunHx = 12 * bsxfun(@minus, trunx, diaga).^2 - 2;
    
    Hx(m, n, n, nobj) = 0;
    for i = 1 : nobj
      Hx(:, i, i, i) = trunHx(:, i);
    end
    
    Hx = Hx + 2 * utils.eyen([m, n, n, nobj], 2, 3);
    
    if m == 1
      Hx = shiftdim(Hx, 1); % (n x n x nobj)
    end
  end
end

