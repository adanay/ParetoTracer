% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [nonlcon] = valnonlcon(nonlcon, doval)
% Validates the nonlinear constraints function.
% The nonlcon function must be either a cell array of function handles or 
% a struct. I.e.,
% - nonlcon = {c, ceq, Jc, Jceq, Hc, Hceq} representing the inequality and 
%   equality constraints together with their respective Jacobians and
%   Hessians. They all can be empty.
% - nonlcon is a struct such that
%   + cx = nonlcon.c(x)
%   + ceqx = nonlcon.ceq(x)
%   + Jcx = nonlcon.Jc(x) 
%   + Jceqx = nonlcon.Jceq(x) 
%   + Hcx = nonlcon.Hc(x)
%   + Hceqx = nonlcon.Hceq(x)

if ~doval 
  % no validation required
  if ~isstruct(nonlcon)
    nonlcon = struct('c', nonlcon{1}, 'ceq', nonlcon{2}, 'Jc', nonlcon{3}, 'Jceq', nonlcon{4}, 'Hc', nonlcon{5}, 'Hceq', nonlcon{6});
  end
  return
end

if iscell(nonlcon)
  l = length(nonlcon);
  
  % constraints
  if l > 0
    c = nonlcon{1};
    if ~isempty(c) && ~ischar(c) && ~isa(c, 'function_handle')
      error('val:valnonlcon:InvalidC', 'The inequality constraints function nonlcon{1} must be a string or a handle or empty.');
    end
  else
    c = [];
  end
  if l > 1
    ceq = nonlcon{2};
    if ~isempty(ceq) && ~ischar(ceq) && ~isa(ceq, 'function_handle')
      error('val:valnonlcon:InvalidCeq', 'The equality constraints function nonlcon{2} must be a string or a handle or empty.');
    end
  else
    ceq = [];
  end
  % Jacobian
  if l > 2
    Jc = nonlcon{3};
    if ~isempty(Jc) && ~ischar(Jc) && ~isa(Jc, 'function_handle')
      error('val:valnonlcon:InvalidJc', 'The inequality Jacobian function nonlcon{3} must be a string or a handle or empty.');
    end
  else
    Jc = [];
  end
  if l > 3
    Jceq = nonlcon{4};
    if ~isempty(Jceq) && ~ischar(Jceq) && ~isa(Jceq, 'function_handle')
      error('val:valnonlcon:InvalidJceq', 'The equality Jacobian function nonlcon{4} must be a string or a handle or empty.');
    end
  else
    Jceq = [];
  end
  %Hessian
  if l > 4
    Hc = nonlcon{5};
    if ~isempty(Hc) && ~ischar(Hc) && ~isa(Hc, 'function_handle')
      error('val:valnonlcon:InvalidHc', 'The inequality Hessian function nonlcon{5} must be a string or a handle or empty.');
    end
  else
    Hc = [];
  end
  if l > 5
    Hceq = nonlcon{6};
    if ~isempty(Hceq) && ~ischar(Hceq) && ~isa(Hceq, 'function_handle')
      error('val:valnonlcon:InvalidHceq', 'The equality Hessian function nonlcon{6} must be a string or a handle or empty.');
    end
  else
    Hceq = [];
  end
elseif isstruct(nonlcon)
  % objective
  if isfield(nonlcon, 'c')
    c = nonlcon.c;
    if ~isempty(c) && ~ischar(c) && ~isa(c, 'function_handle')
      error('val:valnonlcon:InvalidC', 'The inequality constraints function nonlcon.c must be a string or a handle or empty.');
    end
  else
    c = [];
  end
  if isfield(nonlcon, 'ceq')
    ceq = nonlcon.ceq;
    if ~isempty(ceq) && ~ischar(ceq) && ~isa(ceq, 'function_handle')
      error('val:valnonlcon:InvalidCeq', 'The equality constraints function nonlcon.ceq must be a string or a handle or empty.');
    end
  else
    ceq = [];
  end
  % Jacobian
  if isfield(nonlcon, 'Jc')
    Jc = nonlcon.Jc;
    if ~isempty(Jc) && ~ischar(Jc) && ~isa(Jc, 'function_handle')
      error('val:valnonlcon:InvalidJc', 'The inequality Jacobian function nonlcon.Jc must be a string or a handle or empty.');
    end
  else
    Jc = [];
  end
  % Jacobian
  if isfield(nonlcon, 'Jceq')
    Jceq = nonlcon.Jceq;
    if ~isempty(Jceq) && ~ischar(Jceq) && ~isa(Jceq, 'function_handle')
      error('val:valnonlcon:InvalidJceq', 'The equality Jacobian function nonlcon.Jceq must be a string or a handle or empty.');
    end
  else
    Jceq = [];
  end
  %Hessian
  if isfield(nonlcon, 'Hc')
    Hc = nonlcon.Hc;
    if ~isempty(Hc) && ~ischar(Hc) && ~isa(Hc, 'function_handle')
      error('val:valnonlcon:InvalidHc', 'The inequality Hessian function nonlcon.Hc must be a string or a handle or empty.');
    end
  else
    Hc = [];
  end
  if isfield(nonlcon, 'Hceq')
    Hceq = nonlcon.Hceq;
    if ~isempty(Hceq) && ~ischar(Hceq) && ~isa(Hceq, 'function_handle')
      error('val:valnonlcon:InvalidHceq', 'The equality Hessian function nonlcon.Hceq must be a string or a handle or empty.');
    end
  else
    Hceq = [];
  end
else
  c = nonlcon;
  ceq = [];
  Jc = [];
  Jceq = [];
  Hc = [];
  Hceq = [];
  if ~isempty(c) && ~ischar(c) && ~isa(c, 'function_handle')
    error('val:valnonlcon:InvalidNonlcon', 'The function nonlcon must be either a cell array {c, ceq, Jc, Jceq, Hc, Hceq} or a struct with those fields. Alternatively, a function handle c can also be specified.');
  end
end
nonlcon = struct('c', c, 'ceq', ceq, 'Jc', Jc, 'Jceq', Jceq, 'Hc', Hc, 'Hceq', Hceq);
end