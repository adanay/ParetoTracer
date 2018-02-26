% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [vHx, it1, stats] = vheval(v, it0, it1, objfun, lb, ub, multfun, opts, stats)
% Computes the product v' * Hx.
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.

% sizes
nobj = length(it1.fx);

% the Hessian is known to be the identity
if ~isempty(it1.Hident) && it1.Hident
  vHx = v;
  vHx = shiftdim(vHx, -1);
  vHx = repmat(vHx, nobj, 1);
  vHx = permute(vHx, [2 1 3]);
  return
end

% what we have is the Cholesky factors only
if isempty(multfun.vH) && isempty(it1.Hx) && ~isempty(it1.Lx)
  vHx = eval.vlltimes(shiftdim(it1.Lx, -1), v);
  stats.vHCount = stats.vHCount + 1;
  return
end

switch(lower(opts.HessModif))
  case 'off'
    switch(lower(opts.HessApprox))
      case 'fd'
        [vHx, vHCount, ~, ~, Hx, HCount, ~, ~, Hident, ~, fCount, ~, JCount] = eval.vhevalormultorfd(...
         multfun.vH, objfun.H, objfun.f, objfun.J, it1.x, it1.fx, shiftdim(it1.Jx, -1), shiftdim(it1.Hx, -1), v, lb, ub, true, opts.FDForceHess && ~opts.LargeScale, false, opts);
        stats.fCount = stats.fCount + fCount;
        stats.JCount = stats.JCount + JCount;
      case 'off'
        % the Hessian will be approximated by the identity if necessary
        [vHx, vHCount, ~, ~, Hx, HCount, ~, Hident] = eval.vhevalormultorsd(...
         multfun.vH, objfun.H, it1.x, shiftdim(it1.Hx, -1), v, nobj, true, false, false, opts);
      case 'bfgs' 
        [vHx, vHCount, ~, ~, Hx, HCount, ~, ~, Hident] = eval.vhevalormultorqn(...
         multfun.vH, objfun.H, it1.x, shiftdim(it1.Jx, -1), shiftdim(it1.Hx, -1), it0.x, shiftdim(it0.Jx, -1), shiftdim(it0.Hx, -1), v, true, false, false, opts);
    end
    it1.Lx = [];
    it1.Lmodif = [];
  case 'chol'
    switch(lower(opts.HessApprox))
      case 'fd'
        [vHx, vHCount, ~, ~, Hx, HCount, ~, ~, Hident, Lx, Lmodif, ~, fCount, ~, JCount] = eval.vhevalormultorfdandmodchol(...
         multfun.vH, objfun.H, objfun.f, objfun.J, it1.x, it1.fx, shiftdim(it1.Jx, -1), shiftdim(it1.Hx, -1), v, lb, ub, true, opts.FDForceHess && ~opts.LargeScale, false, opts);
        stats.fCount = stats.fCount + fCount;
        stats.JCount = stats.JCount + JCount;
      case 'off'
        % the Hessian will be approximated by the identity if necessary
        [vHx, vHCount, ~, ~, Hx, HCount, ~, Hident, Lx, Lmodif] = eval.vhevalormultorsdandmodchol(...
         multfun.vH, objfun.H, it1.x, shiftdim(it1.Hx, -1), v, nobj, true, false, false, opts);
      case 'bfgs' 
        [vHx, vHCount, ~, ~, Hx, HCount, ~, ~, Hident, Lx, Lmodif] = eval.vhevalormultorqnandmodchol(...
         multfun.vH, objfun.H, it1.x, shiftdim(it1.Jx, -1), shiftdim(it1.Hx, -1), it0.x, shiftdim(it0.Jx, -1), shiftdim(it0.Hx, -1), v, true, false, false, opts);
    end
    it1.Lx = shiftdim(Lx, 1);
    it1.Lmodif = Lmodif;
end

stats.vHCount = stats.vHCount + vHCount;

it1.Hx = shiftdim(Hx, 1);
it1.Hident = Hident;
stats.HCount = stats.HCount + HCount;
end



