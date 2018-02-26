% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [funvals] = valfunvals(funvals, checkW, n, na, naeq, funvalcheck, doval)
% Validates the function values of a point.
% funvals are the known function values of some point x.
% It must be either a cell array or a struct. I.e.,
% - funvals = {fx, Jx, Hx, ax, aeqx, cx, ceqx, dcx, dceqx, Jcx, Jceqx, Hcx, Hceqx} 
% - funvals can also be a structure with those fields.
% They can be empty.
% If checkW is true then,
% - funvals = {fx,  Jx,  Hx,  ax,  aeqx,  cx,  ceqx, dcx, dceqx, Jcx,  Jceqx,  Hcx,  Hceqx,
%          w, wfx, wJx, Hwx, wax, waeqx, wcx, wceqx,            wJcx, wJceqx, Hcwx, Hceqwx} 
% - funvals can also be a structure with those fields.
% - w: Structure with the Lagrange multipliers.
%   + objectives: Objective functions.
%   + lower: Lower bounds.
%   + upper: Upper bounds.
%   + ineqlin: Linear inequalities.
%   + eqlin: Linear equalities.
%   + ineqnonlin: Nonlinear inequalities.
%   + eqnonlin: Nonlinear equalities.

if ~doval 
  % no validation required
  if ~isstruct(funvals)
    if checkW
      wIndex = 14;
      funvals = struct('fx', funvals{1}, 'Jx', funvals{2}, 'Hx', funvals{3}, 'Hident', [], 'Lx', [], 'Lmodif', [],...
                       'ax', funvals{4}, 'aeqx', funvals{5},...
                       'cx', funvals{6}, 'ceqx', funvals{7},...
                       'dcx', funvals{8}, 'dceqx', funvals{9},...
                       'Jcx', funvals{10}, 'Jceqx', funvals{11},...
                       'Hcx', funvals{12}, 'Hceqx', funvals{13},...
                       'w', funvals{14},...
                       'wfx', funvals{wIndex + 1}, 'wJx', funvals{wIndex + 2}, 'Hwx', funvals{wIndex + 3},...
                       'wax', funvals{wIndex + 4}, 'waeqx', funvals{wIndex + 5},...
                       'wcx', funvals{wIndex + 6}, 'wceqx', funvals{wIndex + 7},...
                       'wJcx', funvals{wIndex + 8}, 'wJceqx', funvals{wIndex + 9},...
                       'Hcwx', funvals{wIndex + 10}, 'Hceqwx', funvals{wIndex + 11});
    else
      funvals = struct('fx', funvals{1}, 'Jx', funvals{2}, 'Hx', funvals{3}, 'Hident', [], 'Lx', [], 'Lmodif', [],...
                       'ax', funvals{4}, 'aeqx', funvals{5},...
                       'cx', funvals{6}, 'ceqx', funvals{7},...
                       'dcx', funvals{8}, 'dceqx', funvals{9},....
                       'Jcx', funvals{10}, 'Jceqx', funvals{11},...
                       'Hcx', funvals{12}, 'Hceqx', funvals{13});
    end
  end
  return
end

hasW = false;
if iscell(funvals)
  l = length(funvals);
  
  % objectives
  if l > 0
    fx = funvals{1};
    fx = fx(:)';
  else
    fx = [];
  end
  
  if l > 1
    Jx = funvals{2};
    Jx = Jx(:, :);
  else
    Jx = [];
  end
  
  if l > 2
    Hx = funvals{3};
    Hx = Hx(:, :, :);
  else
    Hx = [];
  end
  Hident = [];
  Lx = [];
  Lmodif = [];
  
  % linear constraints
  if l > 3
    ax = funvals{4};
    ax = ax(:)';
  else
    ax = [];
  end
  
  if l > 4
    aeqx = funvals{5};
    aeqx = aeqx(:)';
  else
    aeqx = [];
  end
  
  % nonlinear constraints
  if l > 5
    cx = funvals{6};
    cx = cx(:)';
  else
    cx = [];
  end
  
  if l > 6
    ceqx = funvals{7};
    ceqx = ceqx(:)';
  else
    ceqx = [];
  end
  
  % square norm inequalities
  if l > 7
    dcx = funvals{8};
  else
    dcx = [];
  end
  
  % square norm equalities
  if l > 8
    dceqx = funvals{9};
  else
    dceqx = [];
  end
  
  if l > 9
    Jcx = funvals{10};
    Jcx = Jcx(:, :);
  else
    Jcx = [];
  end
  
  if l > 10
    Jceqx = funvals{11};
    Jceqx = Jceqx(:, :);
  else
    Jceqx = [];
  end
  
  if l > 11
    Hcx = funvals{12};
    Hcx = Hcx(:, :, :);
  else
    Hcx = [];
  end
  
  if l > 12
    Hceqx = funvals{13};
    Hceqx = Hceqx(:, :, :);
  else
    Hceqx = [];
  end
  
  if checkW
    wIndex = 14;
    
    % Lagrange multipliers
    if l > wIndex - 1
      w = funvals{wIndex};
      if ~isempty(w)
        hasW = true;
      else
        w = struct(...
          'objectives', [],...
          'lower', [],...
          'upper', [],...
          'ineqlin', [],...
          'eqlin', [],...
          'ineqnonlin', [],...
          'eqnonlin', []);
      end
    else
      w = struct(...
        'objectives', [],...
        'lower', [],...
        'upper', [],...
        'ineqlin', [],...
        'eqlin', [],...
        'ineqnonlin', [],...
        'eqnonlin', []);
    end
    
    % weighted objectives
    if hasW && l > wIndex + 0
      wfx = funvals{wIndex + 1};
    else
      wfx = [];
    end
  
    if hasW && l > wIndex + 1
      wJx = funvals{wIndex + 2};
      wJx = wJx(:)';
    else
      wJx = [];
    end
  
    if hasW && l > wIndex + 2
      Hwx = funvals{wIndex + 3};
      Hwx = Hwx(:, :);
    else
      Hwx = [];
    end
  
    % weighted linear constraints
    if hasW && l > wIndex + 3
      wax = funvals{wIndex + 4};
    else
      wax = [];
    end
  
    if hasW && l > wIndex + 4
      waeqx = funvals{wIndex + 5};
    else
      waeqx = [];
    end
  
    % weighted nonlinear constraints
    if hasW && l > wIndex + 5
      wcx = funvals{wIndex + 6};
    else
      wcx = [];
    end
  
    if hasW && l > wIndex + 6
      wceqx = funvals{wIndex + 7};
    else
      wceqx = [];
    end
  
    if hasW && l > wIndex + 7
      wJcx = funvals{wIndex + 8};
      wJcx = wJcx(:)';
    else
      wJcx = [];
    end
  
    if hasW && l > wIndex + 8
      wJceqx = funvals{wIndex + 9};
      wJceqx = wJceqx(:)';
    else
      wJceqx = [];
    end
  
    if hasW && l > wIndex + 9
      Hcwx = funvals{wIndex + 10};
      Hcwx = Hcwx(:, :);
    else
      Hcwx = [];
    end
  
    if hasW && l > wIndex + 10
      Hceqwx = funvals{wIndex + 11};
      Hceqwx = Hceqwx(:, :);
    else
      Hceqwx = [];
    end
  end
elseif isstruct(funvals) % struct
  % objectives
  if isfield(funvals, 'fx')
    fx = funvals.fx;
    fx = fx(:)';
  else
    fx = [];
  end
  
  if isfield(funvals, 'Jx')
    Jx = funvals.Jx;
    Jx = Jx(:, :);
  else
    Jx = [];
  end
  
  if isfield(funvals, 'Hx')
    Hx = funvals.Hx;
    Hx = Hx(:, :, :);
  else
    Hx = [];
  end
  
  if isfield(funvals, 'Hident')
    Hident = funvals.Hident;
  else
    Hident = [];
  end
  
  if isfield(funvals, 'Lx')
    Lx = funvals.Lx;
  else
    Lx = [];
  end
  
  if isfield(funvals, 'Lmodif')
    Lmodif = funvals.Lmodif;
  else
    Lmodif = [];
  end
  
  % linear constraints
  if isfield(funvals, 'ax')
    ax = funvals.ax;
    ax = ax(:)';
  else
    ax = [];
  end
  
  if isfield(funvals, 'aeqx')
    aeqx = funvals.aeqx;
    aeqx = aeqx(:)';
  else
    aeqx = [];
  end
  
  % nonlinear constraints
  if isfield(funvals, 'cx')
    cx = funvals.cx;
    cx = cx(:)';
  else
    cx = [];
  end
  
  if isfield(funvals, 'ceqx')
    ceqx = funvals.ceqx;
    ceqx = ceqx(:)';
  else
    ceqx = [];
  end
  
  % square norm inequalities
  if isfield(funvals, 'dcx')
    dcx = funvals.dcx;
  else
    dcx = [];
  end
  
  % square norm equalities
  if isfield(funvals, 'dceqx')
    dceqx = funvals.dceqx;
  else
    dceqx = [];
  end
  
  if isfield(funvals, 'Jcx')
    Jcx = funvals.Jcx;
    Jcx = Jcx(:, :);
  else
    Jcx = [];
  end
  
  if isfield(funvals, 'Jceqx')
    Jceqx = funvals.Jceqx;
    Jceqx = Jceqx(:, :);
  else
    Jceqx = [];
  end
  
  if isfield(funvals, 'Hcx')
    Hcx = funvals.Hcx;
    Hcx = Hcx(:, :, :);
  else
    Hcx = [];
  end
  
  if isfield(funvals, 'Hceqx')
    Hceqx = funvals.Hceqx;
    Hceqx = Hceqx(:, :, :);
  else
    Hceqx = [];
  end
  
  if checkW
    % Lagrange multipliers
    if isfield(funvals, 'w')
      w = funvals.w;
      hasW = true;
    else
      w = struct(...
        'objectives', [],...
        'lower', [],...
        'upper', [],...
        'ineqlin', [],...
        'eqlin', [],...
        'ineqnonlin', [],...
        'eqnonlin', []);
    end
    
    % weighted objectives
    if hasW && isfield(funvals, 'wfx')
      wfx = funvals.wfx;
    else
      wfx = [];
    end
  
    if hasW && isfield(funvals, 'wJx')
      wJx = funvals.wJx;
      wJx = wJx(:)';
    else
      wJx = [];
    end
  
    if hasW && isfield(funvals, 'Hwx')
      Hwx = funvals.Hwx;
      Hwx = Hwx(:, :);
    else
      Hwx = [];
    end
  
    % weighted linear constraints
    if hasW && isfield(funvals, 'wax')
      wax = funvals.wax;
    else
      wax = [];
    end
  
    if hasW && isfield(funvals, 'waeqx')
      waeqx = funvals.waeqx;
    else
      waeqx = [];
    end
  
    % weighted nonlinear constraints
    if hasW && isfield(funvals, 'wcx')
      wcx = funvals.wcx;
    else
      wcx = [];
    end
  
    if hasW && isfield(funvals, 'wceqx')
      wceqx = funvals.wceqx;
    else
      wceqx = [];
    end
  
    if hasW && isfield(funvals, 'wJcx')
      wJcx = funvals.wJcx;
      wJcx = wJcx(:)';
    else
      wJcx = [];
    end
  
    if hasW && isfield(funvals, 'wJceqx')
      wJceqx = funvals.wJceqx;
      wJceqx = wJceqx(:)';
    else
      wJceqx = [];
    end
  
    if hasW && isfield(funvals, 'Hcwx')
      Hcwx = funvals.Hcwx;
      Hcwx = Hcwx(:, :);
    else
      Hcwx = [];
    end
  
    if hasW && isfield(funvals, 'Hceqwx')
      Hceqwx = funvals.Hceqwx;
      Hceqwx = Hceqwx(:, :);
    else
      Hceqwx = [];
    end
  end
else
  % objectives
  if ~isempty(funvals)
    fx = funvals;
  else
    fx = [];
  end
  Jx = [];
  Hx = [];
  Hident = [];
  Lx = [];
  Lmodif = [];
  % linear constraints
  ax = [];
  aeqx = [];
  % nonlinear constraints
  cx = [];
  ceqx = [];
  dcx = [];  % square norm inequalities
  dceqx = []; % square norm equalities
  Jcx = [];
  Jceqx = [];
  Hcx = [];
  Hceqx = [];
  if checkW
    % Lagrange multipliers
    w = struct(...
        'objectives', [],...
        'lower', [],...
        'upper', [],...
        'ineqlin', [],...
        'eqlin', [],...
        'ineqnonlin', [],...
        'eqnonlin', []);
    % weighted objectives
    wfx = [];
    wJx = [];
    Hwx = [];
    % weighted linear constraints
    wax = [];
    waeqx = [];
    % weighted nonlinear constraints
    wcx = [];
    wceqx = [];
    wJcx = [];
    wJceqx = [];
    Hcwx = [];
    Hceqwx = [];
  end
end

nobj = [];
nc = [];
nceq = [];

% objectives
if ~isempty(fx)
  nobj = length(fx);
  if ~isempty(Jx) && size(Jx, 1) ~= nobj;
    error('val:valfunvals:FxAndJxInconsistentNobj', 'The number of objectives is inconsistent in the specified fx and Jx.');
  end
end
if ~isempty(Jx)
  nobj = size(Jx, 1);
  if size(Jx, 2) ~= n
    error('val:valfunvals:XAndJxInconsistentN', 'The number of variables is inconsistent in the specified x and Jx.');
  end
end
if ~isempty(Hx)
  if size(Hx, 1) ~= n || size(Hx, 2) ~= n
    error('val:valfunvals:XAndHxInconsistentN', 'The number of variables is inconsistent in the specified x and Hx.');
  end
  if ~isempty(nobj) && size(Hx, 3) ~= nobj
    error('val:valfunvals:FxOrJxAndHxInconsistentNObj', 'The number of objectives is inconsistent in the specified fx (or Jx) and Hx.');
  end
end
if ~isempty(Hident) && ~islogical(Hident)
  error('val:valfunvals:NonLogicalHident', 'The Hident field must be a logical value.');
end
if ~isempty(Lx)
  if size(Lx, 1) ~= n || size(Lx, 2) ~= n
    error('val:valfunvals:XAndLxInconsistentN', 'The number of variables is inconsistent in the specified x and Lx.');
  end
  if ~isempty(nobj) && size(Lx, 3) ~= nobj
    error('val:valfunvals:FxOrJxAndLxInconsistentNObj', 'The number of objectives is inconsistent in the specified fx (or Jx) and Lx.');
  end
end
if ~isempty(Lmodif) && ~islogical(Lmodif)
  error('val:valfunvals:NonLogicalLmodif', 'The Lmodif field must be a logical value.');
end

% linear constraints
if ~isempty(ax) && length(ax) ~= na
  error('val:valfunvals:AOrBAndAxInconsistentNa', 'The number of linear inequalities is inconsistent in the specified ax and A (or b).');
end
if ~isempty(aeqx) && length(aeqx) ~= naeq
  error('val:valfunvals:AeqOrBeqAndAeqxInconsistentNaeq', 'The number of linear equalities is inconsistent in the specified aeqx and Aeq (or beq).');
end

% nonlinear inequalities
if ~isempty(cx)
  nc = length(cx);
  if ~isempty(Jcx) && size(Jcx, 1) ~= nc
    error('val:valfunvals:CxAndJcxInconsistentNc', 'The number of objectives is inconsistent in the specified cx and Jcx.');
  end
end
if ~isempty(Jcx)
  nc = size(Jcx, 1);
  if size(Jcx, 2) ~= n
    error('val:valfunvals:XAndJcxInconsistentN', 'The number of variables is inconsistent in the specified x and Jcx.');
  end
end
if ~isempty(Hcx)
  if size(Hcx, 1) ~= n || size(Hcx, 2) ~= n
    error('val:valfunvals:XAndHcxInconsistentN', 'The number of variables is inconsistent in the specified x and Hcx.');
  end
  if ~isempty(nc) && size(Hcx, 3) ~= nc
    error('val:valfunvals:CxOrJcxAndHcxInconsistentNc', 'The number of constraints is inconsistent in the specified cx (or Jcx) and Hcx.');
  end
end

% nonlinear equalities
if ~isempty(ceqx)
  nceq = length(ceqx);
  if ~isempty(Jceqx) && size(Jceqx, 1) ~= nceq
    error('val:valfunvals:CeqxAndJceqxInconsistentNceq', 'The number of objectives is inconsistent in the specified ceqx and Jceqx.');
  end
end
if ~isempty(Jceqx)
  nceq = size(Jceqx, 1);
  if size(Jceqx, 2) ~= n
    error('val:valfunvals:XAndJceqxInconsistentN', 'The number of variables is inconsistent in the specified x and Jceqx.');
  end
end
if ~isempty(Hceqx)
  if size(Hceqx, 1) ~= n || size(Hceqx, 2) ~= n
    error('val:valfunvals:XAndHceqxInconsistentN', 'The number of variables is inconsistent in the specified x and Hceqx.');
  end
  if ~isempty(nceq) && size(Hceqx, 3) ~= nceq
    error('val:valfunvals:CeqxOrJceqxAndHceqxInconsistentNceq', 'The number of constraints is inconsistent in the specified ceqx (or Jceqx) and Hceqx.');
  end
end

% square norm of the constraints
if ~isempty(dcx) && ~isscalar(dcx)
  error('val:valfunvals:NonScalarDcx', 'The square norm of the inequalities dcx must be a scalar.');
end
if ~isempty(dceqx) && ~isscalar(dceqx)
  error('val:valfunvals:NonScalarDceqx', 'The square norm of the equalities dceqx must be a scalar.');
end

if checkW && hasW
  % Lagrange multipliers
  if ~isstruct(w)
    error('val:valfunvals:NonStructW', 'The Lagrange multipliers w must be a struct with the following fields: objectives, lower, upper, ineqlin, eqlin, ineqnonlin, eqnonlin.');
  end
      
  if ~isfield(w, 'objectives')
    w.objectives = [];
  end
  if ~isfield(w, 'lower')
    w.lower = [];
  end
  if ~isfield(w, 'upper')
    w.upper = [];
  end
  if ~isfield(w, 'ineqlin')
    w.ineqlin = [];
  end
  if ~isfield(w, 'eqlin')
    w.eqlin = [];
  end
  if ~isfield(w, 'ineqnonlin')
    w.ineqnonlin = [];
  end
  if ~isfield(w, 'eqnonlin')
    w.eqnonlin = [];
  end
  
  if isempty(w.objectives) && isempty(w.lower) && isempty(w.upper) && isempty(w.ineqlin) && isempty(w.eqlin) && isempty(w.ineqnonlin) && isempty(w.eqnonlin)
    hasW = false;
  else
    % objectives
    if ~isempty(w.objectives)
      if ~isempty(nobj) && length(w.objectives) ~= nobj
        error('val:valfunvals:InconsistentWObjectives', 'The Lagrange multipliers w.objectives are not consistent with the specified number of objectives.');
      end
      if any(w.objectives < 0)
        warning('val:valfunvals:NegativeWObjectives', 'The Lagrange multipliers w.objectives are supposed to be non negative.');
      end
      if ~isempty(wfx) && ~isscalar(wfx)
        error('val:valfunvals:NonScalarWfx', 'The weighted sum of objectives wfx must be a scalar.');
      end
      if ~isempty(wJx) && length(wJx) ~= n
        error('val:valfunvals:XAndWJxInconsistent', 'The number of variables is inconsistent in the specified x and wJx.');
      end
      if ~isempty(Hwx) && (size(Hwx, 1) ~= n || size(Hwx, 2) ~= n)
        error('val:valfunvals:XAndHwxInconsistentN', 'The number of variables is inconsistent in the specified x and Hwx.');
      end
    else
      wfx = [];
      wJx = [];
      Hwx = [];
    end
    
    % lower bounds
    if ~isempty(w.lower)
      if length(w.lower) ~= n
        error('val:valfunvals:InconsistentWLower', 'The Lagrange multipliers w.lower are not consistent with the specified number of variables.');
      end
      if any(w.lower < 0)
        warning('val:valfunvals:NegativeWLower', 'The Lagrange multipliers w.lower are supposed to be non negative.');
      end
    end
    
    % upper bounds
    if ~isempty(w.upper)
      if length(w.upper) ~= n
        error('val:valfunvals:InconsistentWUpper', 'The Lagrange multipliers w.upper are not consistent with the specified number of variables.');
      end
      if any(w.upper < 0)
        warning('val:valfunvals:NegativeWUpper', 'The Lagrange multipliers w.upper are supposed to be non negative.');
      end
    end
    
    % linear inequalities
    if ~isempty(w.ineqlin)
      if length(w.ineqlin) ~= na
        error('val:valfunvals:InconsistentWIneqlin', 'The Lagrange multipliers w.ineqlin are not consistent with the specified number of linear inequalities.');
      end
      if any(w.ineqlin < 0)
        warning('val:valfunvals:NegativeWIneqlin', 'The Lagrange multipliers w.ineqlin are supposed to be non negative.');
      end
      if ~isempty(wax) && ~isscalar(wax)
        error('val:valfunvals:NonScalarWax', 'The weighted sum of linear inequalities wax must be a scalar.');
      end
    else
      wax = [];
    end
    % linear equalities
    if ~isempty(w.eqlin)
      if length(w.eqlin) ~= naeq
        error('val:valfunvals:InconsistentWEqlin', 'The Lagrange multipliers w.eqlin are not consistent with the specified number of linear equalities.');
      end
      if ~isempty(waeqx) && ~isscalar(waeqx)
        error('val:valfunvals:NonScalarWaeqx', 'The weighted sum of linear equalities waeqx must be a scalar.');
      end
    else
      waeqx = [];
    end
    
    % nonlinear inequalities
    if ~isempty(w.ineqnonlin)
      if ~isempty(nc) && length(w.ineqnonlin) ~= nc
        error('val:valfunvals:InconsistentWIneqnonlin', 'The Lagrange multipliers w.ineqnonlin are not consistent with the specified number of nonlinear inequalities.');
      end
      if any(w.ineqnonlin < 0)
        warning('val:valfunvals:NegativeWIneqnonlin', 'The Lagrange multipliers w.ineqnonlin are supposed to be non negative.');
      end
      if ~isempty(wcx) && ~isscalar(wcx)
        error('val:valfunvals:NonScalarWcx', 'The weighted sum of nonlinear inequalities wcx must be a scalar.');
      end
      if ~isempty(wJcx) && length(wJcx) ~= n
        error('val:valfunvals:XAndWJcxInconsistent', 'The number of variables is inconsistent in the specified x and wJcx.');
      end
      if ~isempty(Hcwx) && (size(Hcwx, 1) ~= n || size(Hcwx, 2) ~= n)
        error('val:valfunvals:XAndHcwxInconsistentN', 'The number of variables is inconsistent in the specified x and Hcwx.');
      end
    else
      wcx = [];
      wJcx = [];
      Hcwx = [];
    end
    % nonlinear equalities
    if ~isempty(w.eqnonlin)
      if ~isempty(nceq) && length(w.eqnonlin) ~= nceq
        error('val:valfunvals:InconsistentWEqnonlin', 'The Lagrange multipliers w.eqnonlin are not consistent with the specified number of nonlinear equalities.');
      end
      if ~isempty(wceqx) && ~isscalar(wceqx)
        error('val:valfunvals:NonScalarWceqx', 'The weighted sum of nonlinear equalities wceqx must be a scalar.');
      end
      if ~isempty(wJceqx) && length(wJceqx) ~= n
        error('val:valfunvals:XAndWJceqxInconsistent', 'The number of variables is inconsistent in the specified x and wJceqx.');
      end
      if ~isempty(Hceqwx) && (size(Hceqwx, 1) ~= n || size(Hceqwx, 2) ~= n)
        error('val:valfunvals:XAndHceqwxInconsistentN', 'The number of variables is inconsistent in the specified x and Hceqwx.');
      end
    else
      wceqx = [];
      wJceqx = [];
      Hceqwx = [];
    end
  end
end

if funvalcheck
  val.checkval(fx, 'fx', false);
  if val.checkval(Jx, 'Jx', true)
    Jx = [];
  end
  if val.checkval(Hx, 'Hx', true)
    Hx = [];
  end
  val.checkval(ax, 'ax', false);
  val.checkval(aeqx, 'aeqx', false);
  val.checkval(cx, 'cx', false);
  val.checkval(ceqx, 'ceqx', false);
  val.checkval(dcx, 'dcx', false);
  val.checkval(dceqx, 'dceqx', false);
  if val.checkval(Jcx, 'Jcx', true)
    Jcx = [];
  end
  if val.checkval(Jceqx, 'Jceqx', true)
    Jceqx = [];
  end
  if val.checkval(Hcx, 'Hcx', true)
    Hcx = [];
  end
  if val.checkval(Hceqx, 'Hceqx', true)
    Hceqx = [];
  end
    
  if checkW && hasW
    val.checkval(w.objectives, 'w.objectives', false);
    val.checkval(w.lower, 'w.lower', false);
    val.checkval(w.upper, 'w.upper', false);
    val.checkval(w.ineqlin, 'w.ineqlin', false);
    val.checkval(w.eqlin, 'w.eqlin', false);
    val.checkval(w.ineqnonlin, 'w.ineqnonlin', false);
    val.checkval(w.eqnonlin, 'w.eqnonlin', false);
    
    val.checkval(wfx, 'wfx', false);
    if val.checkval(wJx, 'wJx', true)
      wJx = [];
    end
    if val.checkval(Hwx, 'Hwx', true)
      Hwx = [];
    end
    val.checkval(wax, 'wax', false);
    val.checkval(waeqx, 'waeqx', false);
    val.checkval(wcx, 'wcx', false);
    val.checkval(wceqx, 'wceqx', false);
    if val.checkval(wJcx, 'wJcx', true)
      wJcx = [];
    end
    if val.checkval(wJceqx, 'wJceqx', true)
      wJceqx = [];
    end
    if val.checkval(Hcwx, 'Hcwx', true)
      Hcwx = [];
    end
    if val.checkval(Hceqwx, 'Hceqwx', true)
      Hceqwx = [];
    end
  end
end

if checkW
  funvals = struct('fx', fx, 'Jx', Jx, 'Hx', Hx, 'Hident', Hident, 'Lx', Lx, 'Lmodif', Lmodif,...
                   'ax', ax, 'aeqx', aeqx,...
                   'cx', cx, 'ceqx', ceqx,...
                   'dcx', dcx, 'dceqx', dceqx,...
                   'Jcx', Jcx, 'Jceqx', Jceqx,...
                   'Hcx', Hcx, 'Hceqx', Hceqx,...
                   'w', w,...
                   'wfx', wfx, 'wJx', wJx, 'Hwx', Hwx,...
                   'wax', wax, 'waeqx', waeqx,...
                   'wcx', wcx, 'wceqx', wceqx,...
                   'wJcx', wJcx, 'wJceqx', wJceqx,...
                   'Hcwx', Hcwx, 'Hceqwx', Hceqwx);
else
  funvals = struct('fx', fx, 'Jx', Jx, 'Hx', Hx, 'Hident', Hident, 'Lx', Lx, 'Lmodif', Lmodif,...
                   'ax', ax, 'aeqx', aeqx,...
                   'cx', cx, 'ceqx', ceqx,...
                   'dcx', dcx, 'dceqx', dceqx,...
                   'Jcx', Jcx, 'Jceqx', Jceqx,...
                   'Hcx', Hcx, 'Hceqx', Hceqx);
end
end
