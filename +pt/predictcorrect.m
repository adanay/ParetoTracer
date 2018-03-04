function [it1, itcs, result, stats, stop] = predictcorrect(it0, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, diropts, opts, stats, pcoutfcn)
% One continuation step, i.e., one or more predictors and their respective 
% correctors.
% If diropts = 'IncNegDir', two predictors (in oposite directions) will be 
% computed instead of one for each tangent direction.
% For the bi-objective case, if diropts = 'LeftUp', the tangent direction 
% will be adjusted to minimize the first objective. If nobj > 2, no action 
% will be taken. Analogously, if diropts = 'RightDown', the tangent 
% direction will be adjusted to minimize the second objective.
%
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.
% Assumes that all objective and nonlinear constraint Jacobians it1.Jx,
% it1.Jcx, it1.Jceqx are already computed.
% Assumes that the Lagrange multipliers w are already computed.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% one PC iteration: one or more predictors + their corresponding correctors
stats.PCIts = stats.PCIts + 1;

if ~exist('pcoutfcn', 'var')
  pcoutfcn = [];
end

% stop condition (output function)
if ~isempty(pcoutfcn);
  [stop, ~, it1, ~, ~, stats] = feval(pcoutfcn, 'IT-INIT', it0, it1, [], [], result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, []);   
  if stop
    [stop, ~, it1, ~, ~, stats] = feval(pcoutfcn, 'IT-FINIT', it0, it1, [], [], result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, -1); 
    return
  end
end

% predictor phase
[it1, itps, stats, stop] = pt.predict(it0, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, diropts, opts, stats, pcoutfcn);
if stop
  if ~isempty(pcoutfcn);
    [stop, ~, it1, ~, ~, stats] = feval(pcoutfcn, 'IT-FINIT', it0, it1, itps, [], result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, -1); 
  end
  return;
end

% corrector phase
[it1, itcs, result, stats, stop] = pt.correct(it0, it1, itps, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, pcoutfcn);
if stats.LastOptCount > 0
  stats.LastOptIts = stats.LastOptIts / stats.LastOptCount;
  stats.LastOptDirIts = stats.LastOptDirIts / stats.LastOptCount;
  stats.LastOptLsIts = stats.LastOptLsIts / stats.LastOptCount;
end
if stop
  if ~isempty(pcoutfcn);
    [stop, ~, it1, ~, ~, stats] = feval(pcoutfcn, 'IT-FINIT', it0, it1, itps, itcs, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, -1); 
  end
  return;
end

% stop condition (output function)
if ~isempty(pcoutfcn);
  [stop, ~, it1, ~, ~, stats] = feval(pcoutfcn, 'IT-FINIT', it0, it1, itps, itcs, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, []);   
end
end