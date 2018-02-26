% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [lincon, na, naeq] = vallincon(lincon, n, doval)
% Validates the linear constraints of an optimization solver.
% The lincon must be either a cell array of matrices or a struct. I.e.,
% - lincon = {A, b, Aeq, beq} representing the linear inequality and  
%   equality constraints. They all can be empty.
% - lincon is a struct such that
%   + A = lincon.A
%   + b = lincon.b
%   + Aeq = lincon.Aeq 
%   + beq = lincon.beq

if ~doval 
  % no validation required
  if ~isstruct(lincon)
    lincon = struct('A', lincon{1}, 'b', lincon{2}, 'Aeq', lincon{3}, 'beq', lincon{4});
  end
  na = length(lincon.b);
  naeq = length(lincon.beq);
  return
end 

if isstruct(lincon)
  if(isfield(lincon, 'A'))
    A = lincon.A;
  else
    A = [];
  end
  if(isfield(lincon, 'b'))
    b = lincon.b;
  else
    b = [];
  end
  if(isfield(lincon, 'Aeq'))
    Aeq = lincon.Aeq;
  else
    Aeq = [];
  end
  if(isfield(lincon, 'beq'))
    beq = lincon.beq;
  else
    beq = [];
  end
elseif iscell(lincon)
  l = length(lincon);
  if l > 0
    A = lincon{1};
  else
    A = [];
  end
  if l > 1
    b = lincon{2};
  else
    b = [];
  end
  if l > 2
    Aeq = lincon{3};
  else
    Aeq = [];
  end
  if l > 3
    beq = lincon{4};
  else
    beq = [];
  end
else
  A = lincon;
  b = [];
  Aeq = [];
  beq = [];
end

val.checkval(A, 'A', false);
val.checkval(b, 'b', false);
val.checkval(Aeq, 'Aeq', false);
val.checkval(beq, 'beq', false);

b = b(:);
na = length(b);
if isempty(A)
  A = reshape(A, 0, n);
else
  A = A(:, :);
  [nrow, ncol] = size(A);
  if ncol ~= n
    error('val:vallincon:WrongNumberOfColumnsInA', 'The matrix A must have %d columns. Currently it has %d.', n, ncol);
  end
  if nrow ~= na
    error('val:vallincon:AAndBInconsistent', 'The matrix A must have the same number of rows than b.')
  end
end

beq = beq(:);
naeq = length(beq);
if isempty(Aeq)
  Aeq = reshape(Aeq, 0, n);
else
  Aeq = Aeq(:, :);
  [nrow, ncol] = size(Aeq);
  if ncol ~= n
    error('val:vallincon:WrongNumberOfColumnsInAeq', 'The matrix Aeq must have %d columns. Currently it has %d.', n, ncol);
  end
  if nrow ~= naeq
    error('val:vallincon:AeqAndBeqInconsistent', 'The matrix Aeq must have the same number of rows than beq.')
  end
end

lincon = struct('A', A, 'b', b, 'Aeq', Aeq, 'beq', beq);
end




