% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Lx0, x0, Jx0, x1, Jx1, m, n, nobj] = valqnmodcholargin(Lx0, x0, Jx0, x1, Jx1, opts)
% Validates the arguments of a QN update method.

doval = ~(nargin >= 5 && isfield(opts, 'ValidateInput') && isa(opts.ValidateInput, 'logical') && ~any(opts.ValidateInput));

x0 = x0(:, :);
[m0, n0] = size(x0);

if m0 <= 1
  Jx0 = Jx0(:, :);
else
  Jx0 = Jx0(:, :, :);
end

x1 = x1(:, :);
[m1, n1] = size(x1);

if m1 == 1
  Jx1 = Jx1(:, :);
else
  Jx1 = Jx1(:, :, :);
end

m = max(m0, m1);
n = n1;

% options
opts = val.valqnopts(opts, doval);

if ~doval
  % no validation required
  if m1 == 1
    nobj = size(Jx1, 1);
  else
    nobj = size(Jx1, 2);
  end
  
  if isempty(Lx0)
    if m0 > 1
      Lx0 = utils.eyen([m0 n n nobj], 2, 3); % (m0 x n x n x nobj)
    else
      Lx0 = utils.eyen([n n nobj], 1, 2); % (n x n x nobj)
    end
  end 
  return
end

% validates the iteration points
if isempty(x0) && (~isempty(Jx0) || ~isempty(Lx0))
  error('val:valqnargin:EmptyX0', 'The previous iteration point x0 must be specified since a non empty value of Jx0 and Lx0 was specified.')
end
val.checkval(x0, 'x0', false);
if isempty(x1)
  error('val:valqnargin:EmptyX1', 'The current iteration point x1 must be specified.')
end
val.checkval(x1, 'x1', false);
if m0 > 1 && m ~= 1 && m0 ~= m
  error('val:valqnargin:X0AndX1InconsistentM', 'The number of individuals is inconsistent in the specified x0 and x1.');
end
if n0 ~= n && m0 ~= 0
  error('val:valqnargin:X0AndX1InconsistentN', 'The number of variables is inconsistent in the specified x0 and x1.');
end
if m1 ~= 1 && m ~= 1 && m1 ~= m
  error('val:valqnargin:X0AndX1InconsistentM', 'The number of individuals is inconsistent in the specified x0 and x1.');
end

% validates the Jacobian function values
if isempty(Jx0) && (~isempty(x0) || ~isempty(Lx0))
  error('val:valqnargin:EmptyJx0', 'The previous iteration Jacobian value Jx0 must be specified.')
end
if isempty(Jx1)
  error('val:valqnargin:EmptyJx1', 'The previous iteration Jacobian value Jx1 must be specified.')
end

if m0 <= 1
  [nobjJx0, nJx0] = size(Jx0);
else
  [mJx0, nobjJx0, nJx0] = size(Jx0);
  if mJx0 ~= m0
    error('val:valqnargin:X0AndJx0InconsistentM', 'The number of individuals is inconsistent in the specified x0 and Jx0.');
  end
end
if nJx0 ~= n && m0 ~= 0
  error('val:valqnargin:X0AndJx0InconsistentN', 'The number of variables is inconsistent in the specified Jx0.');
end
if opts.FunValCheck && m0 ~= 0
  val.checkval(Jx0, 'Jx0', false);
end

if m1 == 1
  [nobjJx1, nJx1] = size(Jx1);
else
  [mJx1, nobjJx1, nJx1] = size(Jx1);
  if mJx1 ~= m1
    error('val:valqnargin:X1AndJx1InconsistentM', 'The number of individuals is inconsistent in the specified x1 and Jx1.');
  end
end
if nJx1 ~= n
  error('val:valqnargin:X1AndJx1InconsistentN', 'The number of variables is inconsistent in the specified Jx1.');
end
if opts.FunValCheck
  val.checkval(Jx1, 'Jx1', false);
end

if nobjJx0 ~= nobjJx1 && m0 ~= 0
  error('val:valqnargin:Jx0AndJx1InconsistentNobj', 'The number of objectives is inconsistent in the specified Jx0 and Jx1.');
end

nobj = nobjJx1;

% validates the Hessian function value
if isempty(Lx0)
  if m0 > 1
    Lx0 = utils.eyen([m0 n n nobj], 2, 3); % (m0 x n x n x nobj)
  else
    Lx0 = utils.eyen([n n nobj], 1, 2); % (n x n x nobj)
  end
else
  if m0 == 1
    [n1Lx0, n2Lx0, nobjLx0] = size(Lx0);
  else
    [mLx0, n1Lx0, n2Lx0, nobjLx0] = size(Lx0);
    if mLx0 ~= m0
      error('val:valqnargin:X0AndLx0InconsistentM', 'The number of individuals is inconsistent in the specified x0 and Lx0.');
    end
  end
  if n1Lx0 ~= n || n2Lx0 ~= n
    error('val:valqnargin:X0AndLx0InconsistentN', 'The number of variables is inconsistent in the specified Lx0.');
  end
  if nobjLx0 ~= nobj
    error('val:valqnargin:Lx0InconsistentNobj', 'The number of objectives is inconsistent in the specified Lx0.');
  end
  if opts.FunValCheck
    val.checkval(Lx0, 'Lx0', false);
  end
end
end