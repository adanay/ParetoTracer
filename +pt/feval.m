function [it1, stats] = feval(it1, objfun, lincon, nonlcon, force, opts, stats)
% Ensures that all objectives and constraints of a MOP are evaluated.
% If force is true, the function values will be recomputed even if they are
% not empty.
% Computes also it1.dcx (square norm of the inequalities) and it1.dceqx
% (square norm of the inequalities).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% objectives
if force || isempty(it1.fx)
  [fx, fCount] = eval.feval(objfun.f, it1.x, false, opts);
  it1.fx = fx;
  stats.fCount = stats.fCount + fCount; 
end

% linear inequalities
if force || isempty(it1.ax)
  if ~isempty(lincon.A)
    [ax, aCount] = a(it1.x);
    it1.ax = ax;
    stats.aCount = stats.aCount + aCount; 
  else
    it1.ax = zeros(1, 0);
  end
end
   
% linear equalities
if force || isempty(it1.aeqx)
  if ~isempty(lincon.Aeq)
    [aeqx, aeqCount] = aeq(it1.x);
    it1.aeqx = aeqx;
    stats.aeqCount = stats.aeqCount + aeqCount;
  else
    it1.aeqx = zeros(1, 0);
  end
end
    
% nonlinear inequalities
if force || isempty(it1.cx)
  if ~isempty(nonlcon.c) 
    [cx, cCount] = eval.feval(nonlcon.c, it1.x, true, opts);
    it1.cx = cx;
    stats.cCount = stats.cCount + cCount;
  else
    it1.cx = zeros(1, 0);
  end
end
    
% nonlinear equalities
if force || isempty(it1.ceqx)
  if ~isempty(nonlcon.ceq) 
    [ceqx, ceqCount] = eval.feval(nonlcon.ceq, it1.x, true, opts);
    it1.ceqx = ceqx;
    stats.ceqCount = stats.ceqCount + ceqCount;
  else
    it1.ceqx = zeros(1, 0);
  end
end

% square norm of the inequalities
if force || isempty(it1.dcx)
  it1.dcx = 0;
  if ~isempty(it1.ax)
    it1.dcx = it1.dcx + it1.ax * it1.ax';
  end
  if ~isempty(it1.cx)
    it1.dcx = it1.dcx + it1.cx * it1.cx';
  end
end

% square norm of the equalities
if force || isempty(it1.dceqx)
  it1.dceqx = 0;
  if ~isempty(it1.aeqx)
    it1.dceqx = it1.dceqx + it1.aeqx * it1.aeqx';
  end
  if ~isempty(it1.ceqx)
    it1.dceqx = it1.dceqx + it1.ceqx * it1.ceqx';
  end
end

  function [ax, aCount] = a(x)
    ax = bsxfun(@minus, (lincon.A * x')', lincon.b');
    aCount = size(ax, 1);
  end

  function [aeqx, aeqCount] = aeq(x)
    aeqx = bsxfun(@minus, (lincon.Aeq * x')', lincon.beq');
    aeqCount = size(aeqx, 1);
  end
end

