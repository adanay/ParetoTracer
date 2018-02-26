% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hwvx, it1, stats] = hwveval(w, v, it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, opts, stats)
% Computes the product Hwv = (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v.
% This method does not ensure that Hwx is positive definite (which is 
% supposed to be required in  case that opts.HessModif ~= 'off'). So far 
% only the predictor uses this function and positive definiteness is not a 
% requirement to compute the tangent.

% the Hessian is known to be the identity
if ~isempty(it1.Hident) && it1.Hident
  if ~isempty(w)
    Hwvx = bsxfun(@times, sum(w, 2), v);
  else
    [it1, stats] = pt.kkt(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, opts, stats); % ensures w is computed
    Hwvx = bsxfun(@times, sum(it1.w.objectives, 2), v);
  end
  return
end

% what we have is the Cholesky factors only
if isempty(multfun.Hwv) && (isempty(it1.Hwx) || isempty(w)) && isempty(it1.Hx) && ~isempty(it1.Lx)
  vHx = eval.vlltimes(shiftdim(it1.Lx, -1), v); % (1 x n x nobj)
  stats.vHCount = stats.vHCount + 1;
  
  Hwvx = eval.vhwtimes(vHx, w); % (1 x n) 
  stats.Hwvcount = stats.Hwvcount + 1;
  return
end

if ~isempty(w)
  switch(lower(opts.HessApprox))
    case 'fd'
      [Hwvx, HwvCount, ~, ~, ~, HwCount, ~, ~, vHCount, ~, ~, Hx, HCount, ~, ~, Hident, ~, ~, fx, fCount, ~, wJCount, ~, Jx, JCount] = eval.hwvevalormultorfd(...
        multfun.Hwv, multfun.Hw, multfun.vH, objfun.H, objfun.f, [], objfun.J, it1.x, [], it1.fx, [], shiftdim(it1.Jx, -1), [], [], shiftdim(it1.Hx, -1), w, v, lb, ub, true, opts.FDForceHess && ~opts.LargeScale, opts);
      it1.fx = fx;
      it1.Jx = shiftdim(Jx, 1);
      stats.fCount = stats.fCount + fCount;
      stats.wJCount = stats.wJCount + wJCount;
      stats.JCount = stats.JCount + JCount;
    case 'off'
      % the Hessian will be approximated by the identity if necessary
      [Hwvx, HwvCount, ~, ~, HwCount, ~, ~, vHCount, ~, Hx, HCount, ~, Hident] = eval.hwvevalormultorsd(...
        multfun.Hwv, multfun.Hw, multfun.vH, objfun.H, it1.x, [], [], shiftdim(it1.Hx, -1), w, v, true, false, opts);
    case 'bfgs'
      [Hwvx, HwvCount, ~, ~, HwCount, ~, ~, vHCount, ~, Hx, HCount, ~, Hident] = eval.hwvevalormultorqn(...
        multfun.Hwv, multfun.Hw, multfun.vH, objfun.H, it1.x, shiftdim(it1.Jx, -1), [], [], shiftdim(it1.Hx, -1), it0.x, shiftdim(it0.Jx, -1), shiftdim(it0.Hx, -1), w, v, true, false, opts);
  end
  it1.Hx = shiftdim(Hx, 1);
  it1.Hident = Hident;
  it1.Lx = [];
  it1.Lmodif = [];
  stats.HwvCount = stats.HwvCount + HwvCount;
  stats.HwCount = stats.HwCount + HwCount;
  stats.vHCount = stats.vHCount + vHCount;
  stats.HCount = stats.HCount + HCount;
else % w = it1.w.objectives
  [it1, stats] = pt.kkt(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, opts, stats); % ensures w is computed
  switch(lower(opts.HessApprox))
    case 'fd'
      [Hwvx, HwvCount, ~, ~, Hwx, HwCount, ~, ~, vHCount, ~, ~, Hx, HCount, ~, ~, Hident, wfx, ~, fx, fCount, wJx, wJCount, ~, Jx, JCount] = eval.hwvevalormultorfd(...
        multfun.Hwv, multfun.Hw, multfun.vH, objfun.H, objfun.f, multfun.wJ, objfun.J, it1.x, it1.wfx, it1.fx, it1.wJx, shiftdim(it1.Jx, -1), shiftdim(it1.Hwx, -1), [], shiftdim(it1.Hx, -1), it1.w.objectives, v, lb, ub, true, opts.FDForceHess && ~opts.LargeScale, false, opts);
      it1.wfx = wfx;
      it1.fx = fx;
      it1.wJx = wJx;
      it1.Jx = shiftdim(Jx, 1);
      stats.fCount = stats.fCount + fCount;
      stats.wJCount = stats.wJCount + wJCount;
      stats.JCount = stats.JCount + JCount;
    case 'off'
      % the Hessian will be approximated by the identity if necessary
      [Hwvx, HwvCount, ~, ~, Hwx, HwCount, ~, ~, vHCount, ~, Hx, HCount, ~, Hident] = eval.hwvevalormultorsd(...
        multfun.Hwv, multfun.Hw, multfun.vH, objfun.H, it1.x, shiftdim(it1.Hwx, -1), [], shiftdim(it1.Hx, -1), it1.w.objectives, v, true, false, false, opts);
    case 'bfgs'
      [Hwvx, HwvCount, ~, ~, Hwx, HwCount, ~, ~, vHCount, ~, Hx, HCount, ~, Hident] = eval.hwvevalormultorqn(...
        multfun.Hwv, multfun.Hw, multfun.vH, objfun.H, it1.x, shiftdim(it1.Jx, -1), shiftdim(it1.Hwx, -1), [], shiftdim(it1.Hx, -1), it0.x, shiftdim(it0.Jx, -1), shiftdim(it0.Hx, -1), it1.w.objectives, v, true, false, false, opts);
  end
  it1.Hwx = shiftdim(Hwx, 1);
  it1.Hx = shiftdim(Hx, 1);
  it1.Hident = Hident;
  it1.Lx = [];
  it1.Lmodif = [];
  stats.HwvCount = stats.HwvCount + HwvCount;
  stats.HwCount = stats.HwCount + HwCount;
  stats.vHCount = stats.vHCount + vHCount;
  stats.HCount = stats.HCount + HCount;
end
end

