% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = quad(sizes, doval, a)
% Quadratic problem.

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
    error('mop:misc:quad', 'The vector a of size (n x nobj) should be specified.');
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
  'objectives', nobj);

% objective function, Jacobian, and Hessians
objfun = struct('f', @quadf, 'J', @quadJ, 'H', @quadH);

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
opts.MOPName = 'Quad';

% Pareto set and front
pset = @quadpset;
pfront = [];

  function [fx] = quadf(x)
    x = val.valmopargin(x, [], [], n, nobj, doval);
    
    fx = sum(bsxfun(@minus, x, a).^2, 2);
    fx = permute(fx, [1 3 2]);
  end

  function [Jx] = quadJ(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    Jx = 2 * bsxfun(@minus, x, a);
    Jx = permute(Jx, [1 3 2]);
    
    if m == 1
      Jx = shiftdim(Jx, 1); % (nobj x n)
    end
  end

  function [Hx] = quadH(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    Hx = 2 * utils.eyen([m, n, n, nobj], 2, 3);
         
    if m == 1
      Hx = shiftdim(Hx, 1); % (n x n x nobj)
    end
  end

  function [ps, m] = quadpset(m)
    [S, m] = utils.usimplex(nobj, m); %(m x nobj)
    b = shiftdim(a, 1);
    ps = S * b';
  end
end

