% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hx, fx, fcount, Jx, Jcount, Jundef, e] = hessian(f, J, x, fx, Jx, lb, ub, v, opts)
% Calculates the Hessians of a function using finite differences (FD). This
% function is vectorized, i.e., X may be a matrix where the no of rows
% represents the no of individuals (m) and the no of columns represents the 
% no of variables (n).
%
% f is the objective function.
% J is the Jacobian function.  
% Only one, f or J, can be empty.
%
% x is the evaluating point. For simplicity, if x is a vector, it must be a 
% row vector. Otherwise it will be taken as several individuals of only one 
% variable.
%
% fx is the function value at x. It can be empty. 
% Jx is the Jacobian value at x. It can be empty.
% If more than one individual is specified, fx and Jx must have a size of 
% (m x nobj) and (m x nobj x n), where nobj is the no of objectives and m 
% is the no of individuals. 
% If m = 1, Jx must be of size (nobj x n).
%
% lb and ub specify the valid bounds for changes in variables. They can be
% empty.
%
% v is a direction vector. Set it only if a 2nd order directional  
% derivative is sought. Note that the number of x individuals and v   
% individuals must either coincide or one of them should be one (in order  
% to be repeated), i.e., the no of rows of x and v must either coincide or 
% be one. The no of columns, i.e., the no of variables, must coincide.
%
% opts describes the options for the execution of the algorithm. 
% See fd.defopts.
%
% Returns Hx, the Hessians of f evaluated at x, and fcount and Jcount, the 
% number of evaluations of the function and the Jacobian respectively. Hx
% has a size of (n x n x nobj). If more than one individual is specified, 
% then Hx will have a size of (m x n x n x nobj). 
% If v is nonempty, Hx is the directional derivative in the direction v,
% i.e. v' * Hx (of size 1 x n). In case of multiobjective functions, 
% v' * Hx would have a size of (nobj x n). If more than one individual is 
% specified, then v' * Hx would have a size of (m x nobj x n).
%
% The algorithm also returns fx, Jx, and Jundef. Jundef tells whether a 
% value of the Jacobian is undefined. In this case the algorithm stops and 
% set the Hessian to the identity.    
% e is the step length used to estimate the Hessian.

% default parameters
if nargin < 4
  fx = [];
end
if nargin < 5
  Jx = [];
end
if nargin < 6
  lb = [];
end
if nargin < 7
  ub = [];
end
if nargin < 8
  v = [];
end
if nargin < 9
  opts = [];
end
  
% validation phase
[f, J, x, m, n, fx, Jx, Jundef, lb, ub, v, opts] = val.valfdargin(2, f, J, x, fx, Jx, lb, ub, v, opts);
mx = size(x, 1);

% initialization phase
fcount = 0;
Jcount = 0;
if ~isempty(J)
  if isempty(Jx)
    Jx = jac(x);
  else
    if mx == 1
      Jx = shiftdim(Jx, -1);
    end
  end
  nobj = size(Jx, 2);
else
  if isempty(fx)
    fx = fun(x);
  end
  nobj = size(fx, 2);
end

% step size
et = stepsize();
e = et;

% Hessians approximation phase
if isempty(v)
  if isempty(J)
    hess_f();
  else
    hess_j();
  end
else % directional derivative
  if isempty(J)
    if opts.LargeScale
      hessv_f_lim();
    else
      hessv_f_nonlim();
    end
  else
    hessv_j();
  end
end

if mx == 1 && ~isempty(J) && ~isempty(Jx)
  Jx = shiftdim(Jx, 1);
end
e = et;

  function [] = hess_j()
    % Hessian approximation using Jacobian.
    
    if Jundef
      Hx = eye(n);
      Hx = repmat(shiftdim(Hx, -1), m, 1, 1); % (m x n x n)
      Hx = repmat(Hx, 1, 1, 1, nobj); % (m x n x n x nobj)
    end
    
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
        error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
    end
    
    y = reshape(permute(y, [2 1 3]), n * m, n); % (n * m) x n
    [Jy, Jundef] = jac(y); % (n * m) x nobj x n
    if Jundef
      Hx = eye(n);
      Hx = repmat(shiftdim(Hx, -1), m, 1, 1); % (m x n x n)
      Hx = repmat(Hx, 1, 1, 1, nobj); % (m x n x n x nobj)
    else
      Jy = reshape(Jy, n, m, nobj, n);
    
      e = utils.diag3(e)'; % (n x m)
    
      % do the approximation
      Hx = bsxfun(@rdivide, bsxfun(@minus, Jy, shiftdim(Jx, -1)), e); % (n x m x nobj x n)
      Hx = permute(Hx, [2 4 1 3]); % (m x n x n x nobj)
      Hx = (Hx + permute(Hx, [1 3 2 4])) / 2;
    end
    % if only one individual, the result is of size (n x n x nobj) instead of
    % (1 x n x n x nobj)
    if m == 1
      Hx = shiftdim(Hx, 1);
    end
  end

  function [] = hess_f()
    % Hessian approximation using the objective function.
    
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
        error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
    end
    
    ediag = zeros(n, m); % (n x m)
    fy = zeros(n, m, nobj);
    Hx = zeros(m, n, n, nobj);
    
    for i = 1 : n
      z = y(:, 1 : i, :); % (m x i x n)
      yi = y(:, 1 : i, i); % (m x i)
      ei = e(:, i, i); % (m x 1)
      zi = bsxfun(@plus, yi, ei);
      
      % keeps the point inside the box
      ob = zi(:, i) < lb(i) | zi(:, i) > ub(i); % (m x 1)
      if any(ob(:))
        ei = ~ob .* ei - ob .* ei;
        y(:, i, i) = x(:, 1, i) + ei;
        
        ob = y(:, i, i) < lb(i) | y(:, i, i) > ub(i); % (m x 1)
        if any(ob(:))
          error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
        end

        zi = bsxfun(@plus, yi, ei);
        
        ob = zi(:, i) < lb(i) | zi(:, i) > ub(i); % (m x 1)
        if any(ob(:))
          error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
        end
      end
      
      fy(i, :, :) = fun(reshape(y(:, i, :), m, n));
      
      z(:, :, i) = zi; % (m x i x n)
      z = reshape(permute(z, [2 1 3]), i * m, n); % (i * m) x n
      fz = fun(z); % (i * m) x nobj
      fz = reshape(fz, i, m, nobj); % (i x m x nobj)
      
      e(:, i, i) = ei;
      ediag(i, :) = ei';

      % do the approximation
      H = bsxfun(@rdivide,... 
                 bsxfun(@plus, fz - fy(1 : i, :, :), shiftdim(fx, -1) - fy(i, :, :)),... 
                 bsxfun(@times, ei', ediag(1 : i, :))); % (i x m x nobj)
      Hx(:, 1 : i, i, :) = permute(H, [2 1 3]);
      Hx(:, i, 1 : i, :) = Hx(:, 1 : i, i, :);
    end
    
    % if only one individual, the result is of size (n x n x nobj) instead of
    % (1 x n x n x nobj)
    if m == 1
      Hx = shiftdim(Hx, 1);
    end
  end

  function [] = hessv_j()
    % Second order directional derivative approximation using the Jacobian.
    
    if Jundef
      Hx = repmat(v, 1, 1, nobj);
    end
    
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
        error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
    end
    
    [Jy, Jundef] = jac(y);
    if Jundef
      Hx = repmat(v, 1, 1, nobj);
    else
      % do the approximation
      Hx = bsxfun(@rdivide, bsxfun(@minus, Jy, Jx), e);
    end
    
    % if only one individual, the result is of size (nobj x n) instead of
    % (1 x nobj x n)
    if m == 1
      Hx = shiftdim(Hx, 1);
    end
  end

  function [] = hessv_f_nonlim()
    % Second order directional derivative approximation using the objective 
    % function for non limited memory.
    
    e = bsxfun(@times, e, ones(size(v))); % (m x n)
    e1 = min(e, [], 2); % (m x 1)
    y1 = bsxfun(@plus, x, bsxfun(@times, e1, v)); % (m x n)
    
    % keeps the point inside the box
    ob = any(bsxfun(@lt, y1, lb) | bsxfun(@gt, y1, ub), 2); % (m x 1)
    if any(ob(:))
      e1(ob) = -e1(ob);
      y1 = bsxfun(@plus, x, bsxfun(@times, e1, v));
      
      ob = any(bsxfun(@lt, y1, lb) | bsxfun(@gt, y1, ub), 2);
      if any(ob(:))
        error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
    end
    
    lb = shiftdim(lb, -1); % (1 x 1 x n)
    ub = shiftdim(ub, -1); % (1 x 1 x n)
    
    y1 = permute(y1, [1 3 2]); % (m x 1 x n)
    x = permute(x, [1 3 2]); % (m x 1 x n)
    if mx < m
      x = repmat(x, m, 1, 1);
    end
    e2 = utils.diag3(e); % (m x n x n) 
    y2 = bsxfun(@plus, x, e2); % (m x n x n)
    
    % keeps the point inside the box
    ob = bsxfun(@lt, y2, lb) | bsxfun(@gt, y2, ub); % (m x n x n)
    if any(ob(:))
      e2(ob) = -e2(ob);
      y2 = bsxfun(@plus, x, e2);
      
      ob = bsxfun(@lt, y2, lb) | bsxfun(@gt, y2, ub);
      if any(ob(:))
        error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
    end
    
    z = bsxfun(@plus, y1, e2); % (m x n x n)
    
    % keeps the point inside the box
    ob = bsxfun(@lt, z, lb) | bsxfun(@gt, z, ub); % (m x n x n)
    if any(ob(:))
      e2(ob) = -e2(ob);
      y2 = bsxfun(@plus, x, e2);
      
      ob = bsxfun(@lt, y2, lb) | bsxfun(@gt, y2, ub);
      if any(ob(:))
        error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
      
      z = bsxfun(@plus, y1, e2); % (m x n x n)
      
      ob = bsxfun(@lt, z, lb) | bsxfun(@gt, z, ub); % (m x n x n)
      if any(ob(:))
        error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
    end
    
    fy1 = fun(y1);
    
    y2 = y2(1 : mx, :, :); % (mx x n x n)
    y2 = reshape(permute(y2, [2 1 3]), n * mx, n); % (n * mx) x n
    
    fy2 = fun(y2); % (n * mx) x nobj
    fy2 = reshape(fy2, n, mx, nobj);
    fy2 = bsxfun(@times, fy2, ones(n, m, nobj)); 
    
    z = reshape(permute(z, [2 1 3]), n * m, n); % (n * m) x n
    fz = fun(z); % (n * m) x nobj
    fz = reshape(fz, n, m, nobj); % (n x m x nobj)
    
    e2 = utils.diag3(e2)'; % (n x m) 
    
    % do the approximation
    Hx = bsxfun(@rdivide,...
                bsxfun(@plus, bsxfun(@minus, fz, fy2), shiftdim(bsxfun(@minus, fx, fy1), -1)),...
                bsxfun(@times, e1', e2)); % (n x m x nobj)
    Hx = permute(Hx, [2 3 1]); % (m x nobj x n)
    
    % if only one individual, the result is of size (nobj x n) instead of
    % (1 x nobj x n)
    if m == 1
      Hx = shiftdim(Hx, 1);
    end
  end

  function [] = hessv_f_lim()
    % Second order directional derivative approximation using the objective
    % function for limited memory.
    
    e = bsxfun(@times, e, ones(size(v))); % (m x n)
    e1 = min(e, [], 2); % (m x 1)
    y1 = bsxfun(@plus, x, bsxfun(@times, e1, v)); % (m x n)
    
    % keeps the point inside the box
    ob = any(bsxfun(@lt, y1, lb) | bsxfun(@gt, y1, ub), 2); % (m x 1)
    if any(ob(:))
      e1(ob) = -e1(ob);
      y1 = bsxfun(@plus, x, bsxfun(@times, e1, v));
      
      ob = any(bsxfun(@lt, y1, lb) | bsxfun(@gt, y1, ub), 2);
      if any(ob(:))
        error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
      end
    end
    
    if mx < m
      x = repmat(x, m, 1);
    end
    
    fy1 = fun(y1); % (m x nobj)
    fy2 = zeros(n, m, nobj);
    fz = zeros(n, m, nobj);
    e2 = zeros(n, m);
    
    for i = 1 : n
      y2 = x;
      xi = x(:, i); % (m x 1)
      e2i = e(:, i);
      y2i = xi + e2i;
      
      % keeps the point inside the box
      ob = y2i < lb(i) | y2i > ub(i); % (m x 1)
      if any(ob(:))
        e2i = ~ob .* e2i - ob .* e2i;
        y2i = xi + e2i;
        
        ob = y2i < lb(i) | y2i > ub(i);
        
        if any(ob(:))
          error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
        end
      end
      
      z = y1;
      y1i = y1(:, i); % (m x 1)
      zi = y1i + e2i;
      
      % keeps the point inside the box
      ob = zi < lb(i) | zi > ub(i);
      if any(ob(:))
        e2i = ~ob .* e2i - ob .* e2i;
        y2i = xi + e2i;
        
        ob = y2i < lb(i) | y2i > ub(i);
        if any(ob(:))
          error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
        end
        
        zi = y1i + e2i;
        
        ob = zi < lb(i) | zi > ub(i);
        if any(ob(:))
          error('fd:hessian:InfeasY', 'The box constraints are too restrictive compared to the variation in variables used for FD.')
        end
      end
     
      y2(:, i) = y2i;
      fy2(i, :, :) = bsxfun(@times, fun(y2(1 : mx, :)), ones(m, nobj));
      
      z(:, i) = zi;
      fz(i, :, :) = fun(z);
      
      e2(i, :) = e2i;
    end
    
    % do the approximation
    Hx = bsxfun(@rdivide,...
                bsxfun(@plus, bsxfun(@minus, fz, fy2), shiftdim(bsxfun(@minus, fx, fy1), -1)),...
                bsxfun(@times, e1', e2)); % (n x m x nobj)
    Hx = permute(Hx, [2 3 1]); % (m x nobj x n)
    
    % if only one individual, the result is of size (nobj x n) instead of
    % (1 x nobj x n)
    if m == 1
      Hx = shiftdim(Hx, 1);
    end
  end

  function [fx] = fun(x)
    [fx, tfcount] = eval.feval(f, x, false, opts);
    fcount = fcount + tfcount;
  end

  function [Jx, Jundef] = jac(x)
    [Jx, tJcount, Jundef] = eval.jeval(J, x, true, opts);
    Jcount = Jcount + tJcount;
  end

  function [e] = stepsize()
    isforward = strcmpi(opts.FDType, 'forward');
    
    if isforward && ~isempty(J)
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
