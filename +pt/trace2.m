function [it2, itlast, result, stats, EXITFLAG] = trace2(it0, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, diropts, opts, stats, pcoutfcn)
% Performs a continuation for bi-objective problems in the specified orientation.
% If diropts = 'LeftUp', the algorithm will perform a continuation in the 
% direction to minimize the first objective until it reaches one extreme of
% the Pareto curve. Analogously, if diropts = 'RightDown', the algorithm 
% will perform a continuation in the direction to minimize the second 
% objective until it reaches the other extreme of the Pareto curve.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if ~exist('pcoutfcn', 'var')
  pcoutfcn = [];
end

if strcmpi(diropts, 'LeftUp')
  DIRMIN = 1; % first objective will decrease
  DIRMAX = 2; % second objective will increase
elseif strcmpi(diropts, 'RightDown')
  DIRMIN = 2; % second objective will decrease
  DIRMAX = 1; % first objective will increase
end

it2 = [];
EXITFLAG = [];
while true
  % ensures the KKT multipliers are computed
  [it1, stats] = pt.kkt(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, opts, stats);
  
  % stop condition (PF extreme)
  if abs(1 - it1.w.objectives(DIRMIN)) <= opts.PCEdgeTol
    % last point
    itlast = it1;
    
    EXITFLAG = 1;
    return
  end
  
  % predictor -> corrector
  [it1, itc, result, stats, stop] = pt.predictcorrect(it0, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, diropts, opts, stats, pcoutfcn);
  
  % second point
  if isempty(it2)
    it2 = itc;
  end
  
  % last point
  if isempty(itc)
    itlast = it1;
  else
    itlast = itc;
  end
  
  % stop condition (output function)
  if stop
    EXITFLAG = -1;
    return
  end
  
  % stop condition (no more points found)
  if isempty(itc)
    EXITFLAG = 1;
    return
  end
  
  % stop condition (weak optima)
  if itc.fx(DIRMIN) >= it1.fx(DIRMIN) || itc.fx(DIRMAX) < it1.fx(DIRMAX)
    EXITFLAG = 1;
    return
  end
  
  % stop condition (max number of iterations reached)
  if stats.PCIts == opts.PCMaxIts
    EXITFLAG = 0;
    break
  end
  
  % reset 
  it0 = it1;
  it1 = itc;
end
end

