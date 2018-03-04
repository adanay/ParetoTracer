function [it1, itps, stats, stop] = predict(it0, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, diropts, opts, stats, pcoutfcn)
% Computes the predictor(s) for the Pareto Tracer algorithm.
% It will try to compute one predictor for each tangent direction it1.v.
% If diropts = 'IncNegDir', two predictors (in oposite directions) will be 
% computed instead of one for each tangent direction.
% For the bi-objective case, if diropts = 'LeftUp', the tangent direction 
% will be adjusted to minimize the first objective. If nobj > 2, no action 
% will be taken. Analogously, if diropts = 'RightDown', the tangent 
% direction will be adjusted to minimize the second objective.
% There is no guarantee that there will be exactly one predictor for each
% tangent direction as the algorithm may detect that no progress will be
% done in that direction, e.g., the tangent vector is zero or the predictor 
% is discarded by an output function.
%
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.
% Assumes that all objective and nonlinear constraint Jacobians it1.Jx,
% it1.Jcx, it1.Jceqx are already computed.
% Assumes that the Lagrange multipliers w are already computed.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

stop = false;
discard = false;

%sizes
n = length(it1.x);
nobj = length(it1.fx);

if nobj == 1
  itps = [];
  return;
end

if strcmpi(diropts, 'IncNegDir')
  DIRFLAG = 0;
elseif strcmpi(diropts, 'LeftUp')
  DIRFLAG = 1;
elseif strcmpi(diropts, 'RightDown')
  DIRFLAG = 2;
else
  DIRFLAG = -1;
end

% tangents
[it1, stats] = pt.tangent(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, DIRFLAG, opts, stats);

% number of tangent vectors
m = size(it1.v, 1);
if m == 0
  itps = zeros(0, n);
  return;
end

% no of predictors initially equal to no of tangents (or twice the no of tangents)
if DIRFLAG == 0
  itps(2 * m) = pt.minit();
else
  itps(m) = pt.minit();
end
count = 0; % no of actual predictors

for i = 1 : m 
  v = it1.v(i, :);
 
  % step length
  if it1.vIsSecant
    if isempty(it1.Jvx)
      it1.Jvx = it1.Jx * v';
      stats.JvCount = stats.JvCount + 1;
    end
    normJv = norm(it1.Jvx);
    if normJv == 0
      continue
    end
    it1.t = opts.PCStepObj / normJv; 
  else
    % we already solved an equation in which Jv = d and norm(d) = 1
    it1.t = opts.PCStepObj; 
  end 
  
  v = it1.t * v;
  if norm(v) <= eps
    continue
  end
 
  % predictor
  predictor();
  if stop
    break
  end
  
  % opposite predictor
  if DIRFLAG == 0 
    v = -v;
    predictor();
    if stop
      break
    end
  end
end

itps = itps(1 : count);

  function [] = predictor()
    itp = pt.minit();
    
    % predictor
    itp.x = utils.project(it1.x + v, lb, ub);
    
    % computing the objective value at the predictor
    [fp, fCount] = eval.feval(objfun.f, itp.x, false, opts);
    itp.fx = fp;
    stats.fCount = stats.fCount + fCount; 
    
    % checking the neighborhood
    if ~isempty(result.tree)
      e = opts.PCStepObj * 0.7;
      c = scells.ballcontains(result.tree, result.depth, result.pf, itp.fx, e, opts.LbObj, opts.UbObj);  
      % discard (there is already a corrector in the neighborhood of the predictor)
      if ~isempty(c)
        discard = true;
        return
      end
    end
    
    % stop condition (output function)
    if exist('pcoutfcn', 'var') && ~isempty(pcoutfcn)
      [stop, discard, it1, itp, ~, stats] = pcoutfcn('PRED', it0, it1, itp, [], result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, []);   
      if discard || stop
        return
      end
    end
    
    % QN update of the Hessian at the predictor point such that this info
    % can be passed to the corrector.
    if pt.traceisqn(objfun, multfun, opts)
     
      % computing the Jacobian value at the predictor
      [Jp, JCount, ~, ~, ~, fCount] = eval.jevalorfd(objfun.J, objfun.f, itp.x, itp.fx, lb, ub, true, opts);
      itp.Jx = shiftdim(Jp, 1);
      stats.fCount = stats.fCount + fCount;
      stats.JCount = stats.JCount + JCount; 
     
      % updating the Hessian
      if strcmpi(opts.HessModif, 'off')
        [Hp, ~, ~, ~, Hident] = eval.hevalorqn([], itp.x, shiftdim(itp.Jx, -1), it1.x, shiftdim(it1.Jx, -1), shiftdim(it1.Hx, -1), true, false, opts);
        itp.Hx = shiftdim(Hp, 1);
        itp.Hident = Hident;
        itp.Lx = [];
        itp.Lmodif = [];
      else
        [Hp, ~, ~, ~, Hident, Lp, Lmodif] = eval.hevalorqnandmodchol([], itp.x, shiftdim(itp.Jx, -1), it1.x, shiftdim(it1.Jx, -1), shiftdim(it1.Lx, -1), true, false, opts);
        itp.Hx = shiftdim(Hp, 1);
        itp.Hident = Hident;
        itp.Lx = shiftdim(Lp, 1);
        itp.Lmodif = Lmodif;
      end
    end
    
    count = count + 1;
    itps(count) = itp;
  end
end




