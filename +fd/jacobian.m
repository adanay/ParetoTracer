function [Jx, fx, fcount, e] = jacobian(f, x, fx, lb, ub, v, opts)
% Calculates the Jacobian of a function using finite differences (FD). This
% function is vectorized, i.e., X may be a matrix where the no of rows
% represents the no of individuals (m) and the no of columns represents the 
% no of variables (n).
%
% f is the objective function.
%
% x is the evaluating point. For simplicity, if x is a vector, it must be a 
% row vector. Otherwise it will be taken as several individuals of only one 
% variable.
%
% fx is the function value at x. It can be empty. 
% If more than one individual is specified, fx must have a size of 
% (m x nobj), where nobj is the no of objectives and m is the no of 
% individuals. 
%
% lb and ub specify the valid bounds for changes in variables. They can be
% empty.
%
% v is a direction vector. Set it only if a directional derivative is 
% sought. Note that the number of x individuals and v individuals must  
% either coincide or one of them should be one (in order to be repeated), 
% i.e., the no of rows of x and v must either coincide or be one. The no of
% columns, i.e., the no of variables, must coincide.
%
% opts describes the options for the execution of the algorithm. 
% See fd.defopts.
%
% Returns Jx the Jacobian of f evaluated at x, and fcount the number of
% evaluations of the function. Jx has a size of (nobj x n). If more than 
% one individual is specified, then Jx will have a size of (m x nobj x n). 
% If v is nonempty, Jx is the directional derivative in the direction v,
% i.e., Jx * v (of length nobj). If more than one individual is specified, 
% Jx * v will have a size of (m x nobj).
%
% e is the step length used to estimate the Jacobian.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% default parameters
if nargin < 3
  fx = [];
end
if nargin < 4
  lb = [];
end
if nargin < 5
  ub = [];
end
if nargin < 6
  v = [];
end
if nargin < 7
  opts = [];
end

% validation phase
[f, ~, x, m, n, fx, ~, ~, lb, ub, v, opts] = val.valfdargin(1, f, [], x, fx, [], lb, ub, v, opts);

% initialization phase
fcount = 0;
if isempty(fx)
  fx = fun(x);
end
nobj = size(fx, 2);

% step size
et = stepsize();
e = et;

% Jacobian approximation phase
if isempty(v)
  if opts.LargeScale
    jac_lim();
  else
    jac_nonlim();
  end
else % directional derivative
  jacv();
end

e = et;

  function [] = jac_nonlim()
    % Jacobian approximation for non limited memory.
    
    lb = shiftdim(lb, -1); % (1 x 1 x n)
    ub = shiftdim(ub, -1); % (1 x 1 x n)
    
    x = permute(x, [1 3 2]); % (m x 1 x n)
    e = utils.diag3(e); % (m x n x n)
    y = bsxfun(@plus, x, e); % (m x n x n)
    
    % keeps the point inside the box
    ob = bsxfun(@lt, y, lb) | bsxfun(@gt, y, ub); % (m x n x n)
    if any(ob(:))
      e(ob) = -e(ob);
      y = bsxfun(@plus, x, e);
      
      ob = bsxfun(@lt, y, lb) | bsxfun(@gt, y, ub);
      if any(ob(:))
        error('fd:jacobian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
    end
    
    y = reshape(permute(y, [2 1 3]), n * m, n); % (n * m) x n
    fy = fun(y); % (n * m) x nobj
    fy = reshape(fy, n, m, nobj); 
    
    e = utils.diag3(e)'; % (n x m)
    
    % do the approximation
    Jx = bsxfun(@rdivide, bsxfun(@minus, fy, shiftdim(fx, -1)), e); % (n x m x nobj)
    Jx = permute(Jx, [2 3 1]); % (m x nobj x n)
    
    % if only one individual, the result is of size (nobj x n) instead of
    % (1 x nobj x n)
    if m == 1 
      Jx = shiftdim(Jx, 1);
    end
  end

  function [] = jac_lim()
    % Jacobian approximation for limited memory.
    
    Jx = zeros(m, nobj, n);
    
    for i = 1 : n
      y = x;
      xi = x(:, i);
      ei = e(:, i);
      yi = xi + ei;
      
      % keeps the point inside the box
      ob = yi < lb(i) | yi > ub(i);
      if any(ob(:))
        ei = ~ob .* ei - ob .* ei;
        yi = xi + ei;
        
        ob = yi < lb(i) | yi > ub(i);
        if any(ob(:))
          error('fd:jacobian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
        end
      end
      
      y(:, i) = yi;
      fy = fun(y);
      
      % do the approximation
      Jx(:, :, i) = bsxfun(@rdivide, fy - fx, ei);
    end
    
    % if only one individual, the result is of size (nobj x n) instead of
    % (1 x nobj x n)
    if m == 1
      Jx = shiftdim(Jx, 1);
    end
  end

  function [] = jacv()
    % First order directional derivative approximation.
    
    e = bsxfun(@times, e, ones(size(v))); % (m x n)
    e = min(e, [], 2); % (m x 1)
    y = bsxfun(@plus, x, bsxfun(@times, e, v)); % (m x n)
    
    % keeps the point inside the box
    ob = any(bsxfun(@lt, y, lb) | bsxfun(@gt, y, ub), 2); % (m x 1)
    if any(ob(:))
      e(ob) = -e(ob);
      y = bsxfun(@plus, x, bsxfun(@times, e, v));
      
      ob = any(bsxfun(@lt, y, lb) | bsxfun(@gt, y, ub), 2);
      if any(ob(:))
        error('fd:jacobian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
    end
    
    fy = fun(y);
    
    % do the approximation
    Jx = bsxfun(@rdivide, bsxfun(@minus, fy, fx), e);
    
    % if only one individual, the result is of size (nobj x 1) instead of
    % (1 x nobj)
    if m == 1
      Jx = Jx';
    end
  end

  function [fx] = fun(x)
    [fx, tfcount] = eval.feval(f, x, false, opts);
    fcount = fcount + tfcount;
  end

  function [e] = stepsize()
    isforward = strcmpi(opts.FDType, 'forward');
    
    if isforward
      t = opts.FDStepSize(:)';
    else
      t = opts.FDStepSize2(:)';
    end
    e = bsxfun(@times, t, bsxfun(@max, abs(x), abs(opts.TypicalX(:)')));
    
    if isforward            
      e = e .* ((x >= 0) - (x < 0));                   
    end
    e = sign(e) .* bsxfun(@min, bsxfun(@max, abs(e), opts.FDMinChange(:)'),...
                                opts.FDMaxChange(:)');
  end
end
