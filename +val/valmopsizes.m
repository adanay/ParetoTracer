% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [n, nobj, na, naeq, nc, nceq] = valmopsizes(sizes, doval)
% Validates the sizes passed to a mop builder.
% The sizes must be either a cell array of numbers or 
% a struct. I.e.,
% - sizes = {variables, objectives, ineqlin, eqlin, ineqnonlin, eqnonlin}. 
% - sizes is a struct such that
%   + sizes.variables = n
%   + sizes.objectives = nobj
%   + sizes.ineqlin = na
%   + sizes.eqlin = naeq
%   + sizes.ineqnonlin = nc
%   + sizes.eqnonlin = nceq

if ~doval 
  % no validation required
  if ~isstruct(sizes)
    n = sizes{1};
    nobj = sizes{2};
    na = sizes{3};
    naeq = sizes{4};
    nc = sizes{5};
    nceq = sizes{6};
  else
    n = sizes.variables;
    nobj = sizes.objectives;
    na = sizes.ineqlin;
    naeq = sizes.eqlin;
    nc = sizes.ineqnonlin;
    nceq = sizes.eqnonlin;
  end
  return
end

if iscell(sizes)
  l = length(sizes);
  
  % variables
  if l > 0
    n = sizes{1};
    if ~isa(n, 'double')
      error('val:valmopsizes:InvalidN', 'The number of variables must be a double.');
    end
    if isempty(n)
      n = 0;
    end
  else
    n = 0;
  end
  
  % objectives
  if l > 1
    nobj = sizes{2};
    if ~isa(nobj, 'double')
      error('val:valmopsizes:InvalidNObj', 'The number of objectives must be a double.');
    end
    if isempty(nobj)
      nobj = 0;
    end
  else
    nobj = 0;
  end
  
  % linear inequalities
  if l > 2
    na = sizes{3};
    if ~isa(na, 'double')
      error('val:valmopsizes:InvalidNA', 'The number of linear inequalities must be a double.');
    end
    if isempty(na)
      na = 0;
    end
  else
    na = 0;
  end

  % linear equalities
  if l > 3
    naeq = sizes{4};
    if ~isa(naeq, 'double')
      error('val:valmopsizes:InvalidNAeq', 'The number of linear equalities must be a double.');
    end
    if isempty(naeq)
      naeq = 0;
    end
  else
    naeq = 0;
  end
  
  % nonlinear inequalities
  if l > 4
    nc = sizes{5};
    if ~isa(nc, 'double')
      error('val:valmopsizes:InvalidNC', 'The number of nonlinear inequalities must be a double.');
    end
    if isempty(nc)
      nc = 0;
    end
  else
    nc = 0;
  end

  % nonlinear equalities
  if l > 5
    nceq = sizes{6};
    if ~isa(nceq, 'double')
      error('val:valmopsizes:InvalidNCeq', 'The number of nonlinear equalities must be a double.');
    end
    if isempty(nceq)
      nceq = 0;
    end
  else
    nceq = 0;
  end
elseif isstruct(sizes)
  % variables
  if isfield(sizes, 'variables')
    n = sizes.variables;
    if ~isa(n, 'double')
      error('val:valmopsizes:InvalidN', 'The number of variables must be a double.');
    end
    if isempty(n)
      n = 0;
    end
  else
    n = 0;
  end
 
  % objectives
  if isfield(sizes, 'objectives')
    nobj = sizes.objectives;
    if ~isa(nobj, 'double')
      error('val:valmopsizes:InvalidNObj', 'The number of objectives must be a double.');
    end
    if isempty(nobj)
      nobj = 0;
    end
  else
    nobj = 0;
  end
  
  % linear inequalities
  if isfield(sizes, 'ineqlin')
    na = sizes.ineqlin;
    if ~isa(na, 'double')
      error('val:valmopsizes:InvalidNA', 'The number of linear inequalities must be a double.');
    end
    if isempty(na)
      na = 0;
    end
  else
    na = 0;
  end
  
  % linear equalities
  if isfield(sizes, 'eqlin')
    naeq = sizes.eqlin;
    if ~isa(naeq, 'double')
      error('val:valmopsizes:InvalidNAeq', 'The number of linear equalities must be a double.');
    end
    if isempty(naeq)
      naeq = 0;
    end
  else
    naeq = 0;
  end
  
  % nonlinear inequalities
  if isfield(sizes, 'ineqnonlin')
    nc = sizes.ineqnonlin;
    if ~isa(nc, 'double')
      error('val:valmopsizes:InvalidNC', 'The number of nonlinear inequalities must be a double.');
    end
    if isempty(nc)
      nc = 0;
    end
  else
    nc = 0;
  end
  
  % nonlinear equalities
  if isfield(sizes, 'eqnonlin')
    nceq = sizes.eqnonlin;
    if ~isa(nceq, 'double')
      error('val:valmopsizes:InvalidNCeq', 'The number of nonlinear equalities must be a double.');
    end
    if isempty(nceq)
      nceq = 0;
    end
  else
    nceq = 0;
  end
else
  n = 0;
  nobj = 0;
  na = 0;
  naeq = 0;
  nc = 0;
  nceq = 0;
end
end
