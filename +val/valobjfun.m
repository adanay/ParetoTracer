% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun] = valobjfun(objfun, doval)
% Validates the objective function.
% The objective function must be either a cell array of function handles or 
% a struct. I.e.,
% - objfun = {f, J, H} where f, J, H are function handles that represent the
%   objective, Jacobian, and Hessian functions respectively. They can be
%   empty except f.
% - objfun is a struct such that
%   + fx = objfun.f(x)
%   + Jx = objfun.J(x)
%   + Hx = objfun.H(x)

if ~doval 
  % no validation required
  if ~isstruct(objfun)
    objfun = struct('f', objfun{1}, 'J', objfun{2}, 'H', objfun{3});
  end
  return
end

if iscell(objfun)
  l = length(objfun);
  
  % objective
  if l > 0
    f = objfun{1};
    if isempty(f)
      error('val:valobjfun:EmptyFun', 'An objective function objfun{1} must be specified.');
    end
    if ~ischar(f) && ~isa(f, 'function_handle')
      error('val:valobjfun:InvalidFun', 'The objective function objfun{1} must be a string or a handle.');
    end
  else
    error('val:valobjfun:EmptyFun', 'An objective function objfun{1} must be specified.');
  end
  % Jacobian
  if l > 1
    J = objfun{2};
    if ~isempty(J) && ~ischar(J) && ~isa(J, 'function_handle')
      error('val:valobjfun:InvalidJac', 'The Jacobian function objfun{2} must be a string or a handle or empty.');
    end
  else
    J = [];
  end
  %Hessian
  if l > 2
    H = objfun{3};
    if ~isempty(H) && ~ischar(H) && ~isa(H, 'function_handle')
      error('val:valobjfun:InvalidHess', 'The Hessian function objfun{3} must be a string or a handle or empty.');
    end
  else
    H = [];
  end
elseif isstruct(objfun)
  % objective
  if isfield(objfun, 'f')
    f = objfun.f;
    if isempty(f)
      error('val:valobjfun:EmptyFun', 'An objective function objfun.f must be specified.');
    end
    if ~ischar(f) && ~isa(f, 'function_handle')
      error('val:valobjfun:InvalidFun', 'The objective function objfun.f must be a string or a handle.');
    end
  else
    error('val:valobjfun:EmptyFun', 'An objective function objfun.f must be specified.');
  end
  % Jacobian
  if isfield(objfun, 'J')
    J = objfun.J;
    if ~isempty(J) && ~ischar(J) && ~isa(J, 'function_handle')
      error('val:valobjfun:InvalidJac', 'The Jacobian function objfun.J must be a string or a handle or empty.');
    end
  else
    J = [];
  end
  %Hessian
  if isfield(objfun, 'H')
    H = objfun.H;
    if ~isempty(H) && ~ischar(H) && ~isa(H, 'function_handle')
      error('val:valobjfun:InvalidHess', 'The Hessian function objfun.H must be a string or a handle or empty.');
    end
  else
    H = [];
  end
else
  f = objfun;
  J = [];
  H = [];
  if isempty(f)
    error('val:valobjfun:EmptyObjFun', 'The objective function objfun must be either a cell array {f, J, H} or a struct with those fields. Alternatively, a function handle f can also be specified.');
  end
  if ~ischar(f) && ~isa(f, 'function_handle')
    error('val:valobjfun:InvalidObjFun', 'The objective function objfun must be either a cell array {f, J, H} or a struct with those fields. Alternatively, a function handle f can also be specified.');
  end
end

objfun = struct('f', f, 'J', J, 'H', H);
end