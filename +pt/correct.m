% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [it1, itcs, result, stats, stop] = correct(it0, it1, itps, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, pcoutfcn)
% Computes the corrector(s) for the Pareto Tracer algorithm.
% It will try to compute one corrector for each predictor.
% There is no guarantee that there will be exactly one corrector for each
% predictor as the algorithm may detect that no progress will be done, 
% e.g., there is already a corrector in the neighborhood of that predictor, 
% or the corrector is discarded by some output function.
%
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.
% Assumes that all objective and nonlinear constraint Jacobians it1.Jx,
% it1.Jcx, it1.Jceqx are already computed.
% Assumes that the Lagrange multipliers w are already computed.

stop = false;

% number of predictors
m = length(itps);
if m == 0
  itcs = [];
  return;
end

% no of correctors initially equal to no of predictors
itcs(m) = pt.traceit();
count = 0; % no of actual correctors

for i = 1 : m
  itp = itps(i); 
  
  % corrector
  corrector();
  if stop
    break
  end
end

itcs = itcs(1 : count);

  function [] = corrector()
    itc = pt.traceit();
    
    % corrector
    [minresult, optstats] = pt.minimize(objfun, itp.x, itp, lb, ub, lincon, nonlcon, multfun, opts);
    stats = pt.tracestats(stats, optstats);
    itc = pt.traceit(itc, minresult);
    
    % checking the neighborhood
    if ~isempty(result.tree)
      e = opts.PCStepObj * 0.7; 
      c = scells.ballcontains(result.tree, result.depth, result.pf, itc.fx, e, opts.LbObj, opts.UbObj);      
      % discard (there is already a corrector in the neighborhood)
      if ~isempty(c)
        return
      end
    end
    
    % stop condition (bad multipliers)
    if any(it1.w.objectives >= 1 + opts.PCEdgeTol) || any(it1.w.objectives <= -opts.PCEdgeTol)
      return
    end
    
    % discard (too small or too big steps)
    t = norm(itc.fx - it1.fx) / opts.PCStepObj;
    if ~isempty(it0.x) && (t <= opts.PCMinRelStepObj || t >= opts.PCMaxRelStepObj)
      return
    end
    
    % stop condition (weak optima)
    if ~isempty(it0.x) && (utils.isdominant(it1.fx, itc.fx) || utils.isdominant(itc.fx, it1.fx))
      return
    end
    
    % stop condition (output function)
    if exist('pcoutfcn', 'var') && ~isempty(pcoutfcn)
      [stop, discard, it1, ~, itc, stats] = pcoutfcn('CORRECT', it0, it1, itp, itc, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, []);   
      if discard || stop
        return
      end
    end

    % the corrector is a solution
    [itc, result, stats] = pt.tracesave(itc, opts, result, stats);
    
    % (output function)
    if exist('pcoutfcn', 'var') && ~isempty(pcoutfcn)
      pcoutfcn('SAVE', it0, it1, itp, itc, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, []);   
    end
    
    count = count + 1;
    itcs(count) = itc;
  end
end

