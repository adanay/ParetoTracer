function [it1, stats] = jeval(it1, objfun, lb, ub, nonlcon, force, opts, stats)
% Ensures that all objective and nonlinear constraint Jacobians of a MOP 
% are evaluated.
% If force is true, the Jacobian values will be recomputed even if they are
% not empty.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

n = length(it1.x);

% objectives
if force || isempty(it1.Jx)
  [Jx, JCount, ~, ~, fx, fCount] = eval.jevalorfd(objfun.J, objfun.f, it1.x, it1.fx, lb, ub, true, opts);
  it1.Jx = shiftdim(Jx, 1);
  it1.fx = fx;
  stats.JCount = stats.JCount + JCount;
  stats.fCount = stats.fCount + fCount; 
end
    
% nonlinear inequalities
if force || isempty(it1.Jcx)
  if ~isempty(nonlcon.c) 
    [Jcx, JcCount, ~, ~, cx, cCount] = eval.jevalorfd(nonlcon.Jc, nonlcon.c, it1.x, it1.cx, lb, ub, true, opts);
    it1.Jcx = shiftdim(Jcx, 1);
    it1.cx = cx;
    stats.JcCount = stats.JcCount + JcCount;
    stats.cCount = stats.cCount + cCount; 
  else
    it1.Jcx = zeros(0, n);
  end
end
    
% nonlinear equalities
if force || isempty(it1.Jceqx)
  if ~isempty(nonlcon.ceq) 
    [Jceqx, JceqCount, ~, ~, ceqx, ceqCount] = eval.jevalorfd(nonlcon.Jceq, nonlcon.ceq, it1.x, it1.ceqx, lb, ub, true, opts);
    it1.Jceqx = shiftdim(Jceqx, 1);
    it1.ceqx = ceqx;
    stats.JceqCount = stats.JceqCount + JceqCount;
    stats.ceqCount = stats.ceqCount + ceqCount;
  else
    it1.Jceqx = zeros(0, n);
  end
end
end

