function [f, J, x, m, n, fx, Jx, Jundef, lb, ub, v, opts] = valfdargin(order, f, J, x, fx, Jx, lb, ub, v, opts)
% Validates the arguments of a finite differences approximation method.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

doval = ~(nargin >= 10 && isfield(opts, 'ValidateInput') && isa(opts.ValidateInput, 'logical') && ~any(opts.ValidateInput));

x = x(:, :);
[m, n] = size(x);

if m == 1
  fx = fx(:)';
  Jx = Jx(:, :);
else
  fx = fx(:, :);
  Jx = Jx(:, :, :);
end
Jundef = false;

% box constraints
[lb, ub] = val.valboxcon(lb, ub, n, doval);

% options
opts = val.valfdopts(opts, n, doval);

if ~doval
  % no validation required
  v = v(:, :);
  m = max(m, size(v, 1));
  return
end

% validates the objective and Jacobian functions
if order == 1
  if isempty(f)
    error('val:valfdargin:EmptyFun', 'A function f must be specified.');
  end
elseif order == 2
  if isempty(f) && isempty(J)
    error('val:valfdargin:EmptyFunAndJac', 'Either a function f or a Jacobian J must be specified. Only one can be empty.');
  end
end
if ~isempty(f) && ~ischar(f) && ~isa(f, 'function_handle')
  error('val:valfdargin:InvalidFun', 'The function f must be a string or a handle.');
end
if ~isempty(J) && ~ischar(J) && ~isa(J, 'function_handle')
  error('val:valfdargin:InvalidJac', 'The Jacobian function J must be a string or a handle.');
end

% validates the evaluating point
if isempty(x)
  error('val:valfdargin:EmptyX', 'An evaluating point x must be specified.')
end
val.checkval(x, 'x', false);
ob = bsxfun(@lt, x, lb) | bsxfun(@gt, x, ub);
if any(ob(:))
  error('val:valfdargin:InfeasX', 'The evaluating point x is infeasible.')
end

% validates the vector
if ~isempty(v)
  v = v(:, :);
  [mv, nv] = size(v);
  
  if mv ~= 1 && m ~= 1 && mv ~= m
    error('val:valfdargin:XAndVInconsistentM', 'The number of individuals is inconsistent in the specified x and v.');
  end
  
  if nv ~= n
    error('val:valfdargin:XAndVInconsistentN', 'The number of variables is inconsistent in the specified x and v.');
  end
  val.checkval(v, 'v', false);
end

% validates the objective function value
if ~isempty(fx)
  [mfx, nobj] = size(fx);
  if mfx ~= m
    error('val:valfdargin:XAndFxInconsistent', 'The number of individuals is inconsistent in the specified x and fx.');
  end
  if opts.FunValCheck
    val.checkval(fx, 'fx', false);
  end
else
  nobj = 0;
end

% validates the Jacobian function value
if ~isempty(Jx)
  if m == 1
    [nobjJx, nJx] = size(Jx);
  else
    [mJx, nobjJx, nJx] = size(Jx);
    if mJx ~= m
      error('val:valfdargin:XAndJxInconsistentM', 'The number of individuals is inconsistent in the specified x and Jx.');
    end
  end
  if nJx ~= n
    error('val:valfdargin:XAndJxInconsistentN', 'The number of variables is inconsistent in the specified x and Jx.');
  end
  if ~isempty(fx) && nobjJx ~= nobj
    error('val:valfdargin:FxAndJxInconsistentNobj', 'The number of objectives is inconsistent in the specified fx and Jx.');
  end
  if opts.FunValCheck
    Jundef = val.checkval(Jx, 'Jx', true);
  end
end
if ~isempty(v)
  m = max(m, mv);
end
end