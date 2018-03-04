function [it1, stats] = Mu(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, force, opts, stats)
% Computes the Mu space:
% Hw * Mu = [J' A' Aeq' Jc' Jceq'].
% If force is true, the Mu space will be re-computed even if it is not empty.
% Note that only the active inequalities are considered.
% Also, only the non-active dimensions (respect to the box constraints)
% are considered.
% Mu will have a size of (n x nobj) for unconstrained problems.
% For constrained problems it will have a size of (n x (nobj + N)) 
% where N = na + naeq + nc + nceq is the number of active constraints.
% Here, n refers to the number of non-active dimensions respect to the box 
% constraints. 
%
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.
% Assumes that all objective and nonlinear constraint Jacobians it1.Jx, 
% it1.Jcx, it1.Jceqx are already computed.
% Assumes that the Lagrange multipliers w are already computed.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if ~force && ~isempty(it1.Mu)
  return
end

% ensures active sets are computed
it1 = pt.active(it1, lb, ub, false, opts);

% sizes
n = length(it1.x);
nobj = length(it1.fx);
na = sum(it1.axActive);
naeq = length(it1.aeqx);
nc = sum(it1.cxActive);
nceq = length(it1.ceqx);

% Considering only the active inequality constraints.
Jx = it1.Jx;
A = lincon.A(it1.axActive, :);
Aeq = lincon.Aeq;
Jcx = it1.Jcx(it1.cxActive, :);
Jceqx = it1.Jceqx;
 
% Excluding active box constraints since the tangent will have to be zero  
% in those components.
if it1.r < n % if r == n (corner) we will have to compute the unrestricted tangent and project later
  Jx = Jx(:, ~it1.xActive);
  A = A(:, ~it1.xActive); 
  Aeq = Aeq(:, ~it1.xActive);
  Jcx = Jcx(:, ~it1.xActive);
  Jceqx = Jceqx(:, ~it1.xActive);
end
   
% computing the Mu space        
if pt.tracehasfullhess(objfun, multfun) % full Hessian info
  Mu_nonlim(); 
else
  if pt.tracehashessmult(multfun) % multiply Hessian info 
    Mu_lim(); 
  else % no Hessian info
    if strcmpi(opts.HessApprox, 'off') % no Hessian approximation
      w = 1 +...% sum(it1.w.objectives) +...
          sum(it1.w.ineqlin(it1.axActive)) +...
          sum(it1.w.eqlin) +...
          sum(it1.w.ineqnonlin(it1.cxActive)) +...
          sum(it1.w.eqnonlin);
      if w == 0
        w = 1;
        warning('pt:Mu:BadMultipliers', 'The sum of the KKT multipliers w is 0.');
      end
      it1.Mu = [Jx' A' Aeq' Jcx' Jceqx'] / w;
    elseif strcmpi(opts.HessApprox, 'bfgs') % QN approximation
      Mu_nonlim();
    else % FD approximation
      if opts.LargeScale
        Mu_lim();
      else
        Mu_nonlim();
      end
    end
  end
end

  function [] = Mu_nonlim()
    % Computes the Mu space for non limited memory problems.
    % This will solve a linear system of equations using linsolve: the full 
    % weighted Hessian will be computed.
    
    % Computing the weighted sum of the objective hessians.
    [~, it1, stats] = pt.hweval([], it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, opts, stats);
    
    Hwx = it1.Hwx;
    
    % Excluding active dimensions since the tangent will have to be zero in 
    % those components.
    if it1.r < n % if r == n (corner) we will have to compute the unrestricted tangent and project later
      Hwx = it1.Hwx(~it1.xActive, ~it1.xActive);
    end
    
    % linear system       
    o.SYM = true;
    [Mu, R] = linsolve(Hwx, [Jx' A' Aeq' Jcx' Jceqx'], o);
    if R == 0 || val.checkval(Mu, 'Mu', true)
      w = 1 +...% sum(it1.w.objectives) +...
          sum(it1.w.ineqlin(it1.axActive)) +...
          sum(it1.w.eqlin) +...
          sum(it1.w.ineqnonlin(it1.cxActive)) +...
          sum(it1.w.eqnonlin);
      if w == 0
        w = 1;
        warning('pt:Mu:BadMultipliers', 'The sum of the KKT multipliers w is 0.');
      end
      it1.Mu = [Jx' A' Aeq' Jcx' Jceqx'] / w;
    else
      it1.Mu = Mu;
    end
    it1.HwxRcond = R;
%     if R < nobj % this is not correct: R may be the condition number 
%       warning('pt:Mu:rankDeficientMatrix', 'The weighted sum of Hessians Hw is rank deficient with rank=%d and nobj=%d.', R, nobj);
%     end
  end

  function [] = Mu_lim()
    % Computes the Mu space for limited memory problems.
    % Instead of computing the full weighted Hessian, a vector Hwvx will be 
    % computed at each iteration.
    
    Mu(n, nobj + na + naeq + nc + nceq) = 0;
    
    offset = 0;
    npcg(Jx, nobj, offset);
    
    offset = offset + nobj;
    npcg(A, na, offset);
    
    offset = offset + na;
    npcg(Aeq, naeq, offset);
    
    offset = offset + naeq;
    npcg(Jcx, nc, offset);
    
    offset = offset + nc;
    npcg(Jceqx, nceq, offset);

    it1.Mu = Mu;
    
    function [] = npcg(B, nb, offset)
      for i = 1 : nb
        [Mu(:, offset + i), SUBEXITFLAG, ~, ~, ~] = pcg(@subHwv, B(i, :)', [], opts.PCGMaxIts);
        if SUBEXITFLAG ~= 0
          warning('pt:Mu:MuSubProbBadTermination', 'The optimal solution for the Mu space was not found.');
        end
      end
    end
  end

  function [Hwvx] = subHwv(v)
    % Multiply function to solve the linear system of equations 
    % Hw * Mu = [J' A' Aeq' Jc' Jceq'].
    
    [Hwvx, it1, stats] = pt.hwveval([], v(:)', it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, opts, stats);
    
    % Excluding active constraints since the tangent will have to be zero in 
    % those components.
    if it1.r < n % if r == n (corner) we will have to compute the unrestricted tangent and project later
      Hwvx = Hwvx(:, ~it1.xActive);
    end
    
    Hwvx = Hwvx';
  end
end

