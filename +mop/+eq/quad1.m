% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = quad1(sizes, doval, a)
% Quadratic problem with one linear equality constraint.

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

sizes = struct(...
  'variables', n,... 
  'objectives', nobj,...
  'ineqlin', 0,...
  'eqlin', 1,...
  'ineqnonlin', 0,...
  'eqnonlin', 0);

% objective function, Jacobian, and Hessians
objfun = struct('f', @quad1f, 'J', @quad1J, 'H', @quad1H);

% box constraints
lb = [];
ub = [];

% linear constraints
lincon = struct(...
'A', [],...
'b', [],...
'Aeq', [0.5, -1, zeros(1, n - 2)],...
'beq', 0);

% nonlinear constraints
nonlcon = [];      
           
% multiply functions
multfun = [];

% options
opts = pt.defopts(n, nobj);
opts.UseVectorized = true;
opts.MOPName = 'Quad prob with lin eq';

% Pareto set and front
pset = [];
pfront = [];

  function [fx] = quad1f(x)
    x = val.valmopargin(x, [], [], n, nobj, doval);
    
    fx = sum(bsxfun(@minus, x, a).^2, 2);
    fx = permute(fx, [1 3 2]);
  end

  function [Jx] = quad1J(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    Jx = 2 * bsxfun(@minus, x, a);
    Jx = permute(Jx, [1 3 2]);
    
    if m == 1
      Jx = shiftdim(Jx, 1); % (nobj x n)
    end
  end

  function [Hx] = quad1H(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    Hx = 2 * utils.eyen([m, n, n, nobj], 2, 3);
         
    if m == 1
      Hx = shiftdim(Hx, 1); % (n x n x nobj)
    end
  end
end

