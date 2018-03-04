function [it1, stats, dirargout] = kkt(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, force, opts, stats, kktopts)
% Ensures the KKT multipliers are computed.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if ~all(it1.w.objectives == 0) && ~force
  return
end

if nargin < 13 || isempty(kktopts)
  kktopts = optimoptions(@lsqnonlin,...
    'Algorithm', 'levenberg-marquardt',...
    'MaxIter', opts.OptDirMaxIts,...
    'Display', 'off');
end

[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(it1);

% ensures active sets are computed
it1 = pt.active(it1, lb, ub, false, opts);
rLb = it1.rLb;
rUb = it1.rUb;

N = nobj + na + naeq + nc + nceq + rLb + rUb;

% wLb = zeros(1, N);
% if naeq > 0
%   wLb(nobj + na : nobj + na + naeq) = -Inf;
% end
% if nceq > 0
%   wLb(nobj + na + naeq + nc : nobj + na + naeq + nc + nceq) = -Inf;
% end
% wUb = Inf(1, N);

wLb = [];
wUb = [];
  
[w, resnorm, residual, EXITFLAG, output, lambda] = lsqnonlin(@KKT, zeros(1, N), wLb, wUb, kktopts);
dirargout = {w, resnorm, residual, EXITFLAG, output, lambda};

% Lagrange multipliers
it1.w.objectives = w(1 : nobj);
i = nobj;
it1.w.ineqlin = w(i + 1 : i + na);
i = i + na;
it1.w.eqlin = w(i + 1 : i + naeq);
i = i + naeq;
it1.w.ineqnonlin = w(i + 1 : i + nc);
i = i + nc;
it1.w.eqnonlin = w(i + 1 : i + nceq);
i = i + nceq;
it1.w.lower = zeros(1, n);
it1.w.lower(it1.xLbActive) = w(i + 1 : i + rLb);
i = i + rLb;
it1.w.upper = zeros(1, n);
it1.w.upper(it1.xUbActive) = w(i + 1 : i + rUb);

  function [Fw] = KKT(w)
    Fw = [JLagrange(w(:)'),...
          sum(w(1 : nobj)) - 1];
  end

  function [JLx] = JLagrange(w)
    i = 0;
    
    % objectives
    [wJx, wJCount, ~, ~, Jx, JCount, ~, ~, ~, ~, fx, fCount] = eval.wjevalormultorfd([], objfun.J, objfun.f, it1.x, [], it1.fx, shiftdim(it1.Jx, -1), w(i + 1 : i + nobj), lb, ub, true, false, opts); 
    it1.fx = fx;
    it1.Jx = shiftdim(Jx, 1);
    stats.fCount = stats.fCount + fCount;
    stats.wJCount = stats.wJCount + wJCount;
    stats.JCount = stats.JCount + JCount;
    i = i + nobj;
    
    % linear inequalities
    if na == 0
      wJax = zeros(1, n);
    else
      wJax = sum(bsxfun(@times, w(i + 1 : i + na)', lincon.A));
    end
    i = i + na;
    
    % linear equalities
    if naeq == 0
      wJaeqx = zeros(1, n);
    else
      wJaeqx = sum(bsxfun(@times, w(i + 1 : i + naeq)', lincon.Aeq));
    end
    i = i + naeq;
    
    % nonlinear inequalities
    if nc == 0
      wJcx = zeros(1, n);
    else
      [wJcx, wJcCount, ~, ~, Jcx, JcCount, ~, ~, ~, ~, cx, cCount] = eval.wjevalormultorfd([], nonlcon.Jc, nonlcon.c, it1.x, [], it1.cx, shiftdim(it1.Jcx, -1), w(i + 1 : i + nc), lb, ub, true, false, opts); 
      it1.cx = cx;
      it1.Jcx = shiftdim(Jcx, 1);
      stats.cCount = stats.cCount + cCount;
      stats.wJcCount = stats.wJcCount + wJcCount;
      stats.JcCount = stats.JcCount + JcCount; 
    end
    i = i + nc;
    
    % nonlinear equalitites
    if nceq == 0
      wJceqx = zeros(1, n);
    else
      [wJceqx, wJceqCount, ~, ~, Jceqx, JceqCount, ~, ~, ~, ~, ceqx, ceqCount] = eval.wjevalormultorfd([], nonlcon.Jceq, nonlcon.ceq, it1.x, it1.wceqx, it1.ceqx, shiftdim(it1.Jceqx, -1), w(i + 1 : i + nceq), lb, ub, true, false, opts); 
      it1.ceqx = ceqx;
      it1.Jceqx = shiftdim(Jceqx, 1);
      stats.ceqCount = stats.ceqCount + ceqCount;
      stats.wJceqCount = stats.wJceqCount + wJceqCount;
      stats.JceqCount = stats.JceqCount + JceqCount;
    end
    i = i + nceq;
    
    % active box constraints
    wI = zeros(1, n);
    if rLb > 0
      wI(it1.xLbActive) = w(i + 1 : i + rLb);  
    end
    i = i + rLb;
    if rUb > 0
      wI(it1.xUbActive) = w(i + 1 : i + rUb);  
    end
    
    JLx = wJx + wJax + wJaeqx + wJcx + wJceqx + wI; % (1 x n)
  end
end

