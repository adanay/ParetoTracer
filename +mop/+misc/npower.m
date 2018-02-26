% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = npower(sizes, doval, a, gamma)
% Norm to the power of gamma.

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
    error('mop:misc:npower', 'The vector a of size (n x nobj) should be specified.');
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
  gamma = 0.5;
end

sizes = struct(...
  'variables', n,...
  'objectives', nobj);

% objective function, Jacobian, and Hessians
objfun = struct('f', @npowerf, 'J', @npowerJ, 'H', @npowerH);

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
opts.MOPName = 'Power';

% Pareto set and front
pset = @npowerpset;
pfront = [];

  function [fx] = npowerf(x)
    x = val.valmopargin(x, [], [], n, nobj, doval);
    
    X = permute(bsxfun(@minus, x, a), [1 2 4 3]);
    fx = utils.normn(X, 2, 3).^gamma; % (m x 1 x 1 x nobj)
    fx = permute(fx, [1 4 3 2]);
  end

  function [Jx] = npowerJ(x)
    [x, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    X = permute(bsxfun(@minus, x, a), [1 2 4 3]);
    
    q1 = gamma * utils.normn(X, 2, 3).^(gamma - 2); % (m x 1 x 1 x nobj)
    q1 = permute(q1, [1 2 4 3]); % (m x 1 x nobj)
    q2 = bsxfun(@minus, x, a); % (m x n x nobj)
    
    Jx = bsxfun(@times, q1, q2); % (m x n x nobj)
    Jx = permute(Jx, [1 3 2]);
    
    if m == 1
      Jx = shiftdim(Jx, 1); % (nobj x n)
    end
  end

  function [Hx] = npowerH(x)
    [~, ~, ~, m] = val.valmopargin(x, [], [], n, nobj, doval);
    
    X = permute(bsxfun(@minus, x, a), [1 2 4 3]);
    N = utils.normn(X, 2, 3);
    
    q1 = gamma * (gamma - 2) * N.^(gamma - 4); % (m x 1 x 1 x nobj)
    q1 = permute(q1, [1 2 4 3]); % (m x 1 x nobj)
    q1 = permute(q1, [1 2 4 3]);
    
    q2 = bsxfun(@minus, x, a); % (m x n x nobj)
    q2 = permute(q2, [1 2 4 3]); % (m x n x 1 x nobj)
    q2 = qqtimes(q2);
    
    s1 = bsxfun(@times, q1, q2); % (m x n x n x nobj)
    
    q1 = gamma * N.^(gamma - 2); % (m x 1 x 1 x nobj)
    q2 = utils.eyen([m, n, n, nobj], 2, 3);
    s2 = bsxfun(@times, q1, q2); % (m x n x n x nobj)
    
    Hx = s1 + s2;
    
    if m == 1
      Hx = shiftdim(Hx, 1); % (n x n x nobj)
    end
  end

  function [qq] = qqtimes(q)
    m = size(q, 1);
    
    qq(m, n, n, nobj) = 0;
    for i = 1 : m
      for j = 1 : nobj
        qij = q(i, :, :, j); % (1 x n)
        qqij = qij' * qij; % (n x n)
        qq(i, :, :, j) = qqij;
      end
    end
  end

  function [ps, m] = npowerpset(m)
    [S, m] = utils.usimplex(nobj, m); %(m x nobj)
    b = shiftdim(a, 1);
    ps = S * b';
  end
end