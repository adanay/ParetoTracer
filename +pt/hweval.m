% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hwx, it1, stats] = hweval(w, it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, force, opts, stats)
% Computes the product Hw = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj.
% If force is true, the weighted Hessian will be recomputed even if it is
% not empty.
% This method does not ensure that Hwx is positive definite (which is 
% supposed to be required in  case that opts.HessModif ~= 'off'). So far 
% only the predictor uses this function and positive definiteness is not a 
% requirement to compute the tangent.

if isempty(w) && ~isempty(force) && ~force && ~isempty(it1.Hwx)
  return
end

if isempty(w)
  [it1, stats] = pt.kkt(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, opts, stats);
  w = it1.w.objectives;
  persists = true;
else
  persists = false;
end

% the Hessian is known to be the identity
if ~isempty(it1.Hident) && it1.Hident
  Hwx = bsxfun(@times, sum(w, 2), shiftdim(eye(length(it1.x)), -1));
  if persists
    it1.Hwx = shiftdim(Hwx, 1);
  end
  return
end

% what we have is the Cholesky factors only
if isempty(multfun.Hw) && isempty(it1.Hx) && ~isempty(it1.Lx)
  Hx = eval.lltimes(shiftdim(it1.Lx, -1));
  it1.Hx = shiftdim(Hx, 1);
  it1.Lmodif = false;
end

switch(lower(opts.HessApprox))
  case 'fd'
    [Hwx, HwCount, ~, ~, Hx, HCount, ~, ~, Hident, wfx, ~, fx, fCount, wJx, wJCount, ~, Jx, JCount] = eval.hwevalormultorfd(...
      multfun.Hw, objfun.H, objfun.f, [], objfun.J, it1.x, it1.wfx, it1.fx, it1.wJx, shiftdim(it1.Jx, -1), shiftdim(it1.Hx, -1), w, lb, ub, true, opts.FDForceHess && ~opts.LargeScale, false, opts);
    if persists
      it1.wfx = wfx;
      it1.wJx = wJx;
    end
    it1.fx = fx;
    it1.Jx = shiftdim(Jx, 1);
    stats.fCount = stats.fCount + fCount;
    stats.wJCount = stats.wJCount + wJCount;
    stats.JCount = stats.JCount + JCount;
  case 'off'
    [Hwx, HwCount, ~, ~, Hx, HCount, ~, Hident] = eval.hwevalormultorsd(multfun.Hw, objfun.H, it1.x, shiftdim(it1.Hx, -1), w, true, false, false, opts);
  case 'bfgs'
    [Hwx, HwCount, ~, ~, Hx, HCount, ~, ~, Hident] = eval.hwevalormultorqn(multfun.Hw, objfun.H, it1.x, shiftdim(it1.Jx, -1), shiftdim(it1.Hx, -1), it0.x, shiftdim(it0.Jx, -1), shiftdim(it0.Hx, -1), w, true, false, false, opts);
end

if persists
  it1.Hwx = shiftdim(Hwx, 1);
end
stats.HwCount = stats.HwCount + HwCount;

it1.Hx = shiftdim(Hx, 1);
it1.Hident = Hident;
it1.Lx = [];
it1.Lmodif = [];
stats.HCount = stats.HCount + HCount; 
end

