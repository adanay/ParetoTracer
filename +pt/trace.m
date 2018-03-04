function [result, stats, EXITFLAG] = trace(objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts)
% Pareto Tracer (PT). A continuation method for multiobjective optimization
% problems (MOPs). Pareto Tracer performs a continuation on the set of (local)
% optima of a constrained MOP starting at a given (solution) point. 
%
% This version supports only equality constraints and box constrains.
% Support for general inequality constraints is still under construction.
% If some inequality constraint is passed, it is ignored.
%
% pt.trace attempts to trace the set of (local) optima (the Pareto set and 
% the Pareto front) of problems of the form:
%    min f(x)    subject to: A*x <= b, Aeq*x  = beq (linear constraints)
%     x                     c(x) <= 0, ceq(x) = 0   (nonlinear constraints)
%                             lb <= x <= ub         (box contraints)
% pt.trace implements a predictor corrector (PC) approach in which at each
% iteration the tangent space to the Pareto set is computed plus a set of  
% predictor points distributed around the current solution. By utilizing  
% pt.minimize, those predictors are corrected to find new solutions. The  
% process is repeated at each new solution until a stopping criteria is  
% satisfied.
%
% The input parameters are:
%
% objfun: It must be either a cell array of function handles or a struct. 
% I.e., objfun = {f, J, H} where f, J, H are function handles that represent 
% the objective, Jacobian, and Hessian functions respectively. They can be
% empty except f. J and H will be approximated if not provided.
% - objfun can also be a function handle. In this case, it will be assumed 
%   to be the objective function only. I.e., objfun = f. 
% - f is a function handle of the form y = f(x) where x is a vector of n 
%   components and y is a vector of nobj components.
% - J is a function handle of the form y = J(x) where x is a vector of n 
%   components and y is a matrix of size (nobj x n).
% - H is a function handle of the form y = H(x) where x is a vector of n 
%   components and y is a block matrix of size (n x n x nobj).
%
% x0: The initial guess. It must be a vector of n dimensions.
%
% funvals0: The known function values at x0.
% It must be either a cell array or a struct. I.e.,
% - funvals0 = {fx, Jx, Hx, ax, aeqx, cx, ceqx, dcx, dceqx, Jcx, Jceqx} 
% - funvals0 can also be a structure with those fields.
% - fx = f(x0), Jx = J(x0), Hx = H(x0).
% - ax = A * x0 - b, aeqx = Aeq * x0 - beq.
% - cx = c(x0), ceqx = ceq(x0).
% - dcx = norm(cx)^2, dceqx = norm(ceqx)^2.
% - Jcx = Jc(x0), Jceqx = Jceq(x0).
% They all can be empty.
%
% lb and ub: Vectors that represent the box constraints. They must have n
% components.
%
% lincon: It must be either a cell array of matrices or a struct. I.e.,
% - lincon = {A, b, Aeq, beq} representing the linear inequality and  
%   equality constraints. 
% - lincon can also be a structure with those fields.
% They all can be empty.
%
% nonlcon: It must be either a cell array of function handles or a struct. I.e.,
% - nonlcon = {c, ceq, Jc, Jceq} representing the inequality and 
%   equality constraints together with their respective Jacobians. 
% - nonlcon can also be a structure with those fields.
% They all can be empty. If the Jacobians are not provided, they will be
% approximated.
%
% multfun: It must be either a cell array of function handles or a struct. I.e.,
% - multfun = {vH, Hw, Hwv} representing the Hessian multiply functions. 
% - multfun can also be a structure with those fields.
% - vH is a function handle of the form y = vH(x, v) where
%   y = [v' * H1; v' * H2; ...; v' * Hnobj].
%   The result y has a size of (nobj x n).
% - Hw is a function handle of the form y = Hw(x, w) where 
%   y = H1 * w1 + H2 * w2 + ... + Hnobj * wnobj. 
%   The result y has a size of (n x n), i.e., the weighted sum of Hessians.
% - Hwv is a function handle of the form y = Hwv(x, w, v) where 
%   y = (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v.
%   The result y is a vector of n components.
% They all can be empty.
%
% opts: Options for the execution of the algorithm. See pt.defopts.
%
% The output parameters are:
%
% result: Describes the solution. It is a struct with the following fields:
% - ps: The (local) Pareto set of the problem: A matrix of size (m x n) where 
%   m is the number of solutions found and n is the number of variables.
% - pf: The (local) Pareto front of the problem: A matrix of size (m x nobj)  
%   where m is the number of solutions found and nobj is the number of objectives.
%
% EXITFLAG: Describes the exit condition.
%  0: Number of iterations exceeded opts.PCMaxIts.
%  1: No more solution points found.
% -1: Stopped by an output function.
%
% stats: Statistics.
% - Count: Number of solutions found.
% - PCIts: Number of iterations taken by the continuation algorithm, where 
%   each iteration consists of 
%   + one predictor stage (one or more predictor vectors are computed),
%   + and one corrector (minimization) stage (those predictors are corrected).
% - OptIts: Average iterations taken by the corrector (minimization) 
%   algorithm.
% - OptDirIts: Average iterations taken by the corrector direction 
%   subproblem.
% - OptLsIts: Average backtrack iterations taken by the corrector line 
%   search to satisfy the Armijo condition.
%
% - fCount: Number of function evaluations.
% - JCount: Number of Jacobian evaluations.
% - HCount: Number of Hessian evaluations.
% - vHCount, HwCount, HwvCount: Multiply function evaluations.
% - aCount: Number of linear inequality constraint evaluations (A * x - b).
% - aeqCount: Number of linear equality constraint evaluations (Aeq * x - beq).
% - cCount: Number of nonlinear inequality constraint evaluations.
% - ceqCount: Number of nonlinear equality constraint evaluations.
% - JcCount: Number of nonlinear inequality constraint Jacobian evaluations.
% - JceqCount: Number of nonlinear equality constraint Jacobian evaluations.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% default parameters
if nargin < 4
  lb = [];
end
if nargin < 5
  ub = [];
end
if nargin < 6
  lincon = [];
end
if nargin < 7
  nonlcon = [];
end
if nargin < 8
  multfun = [];
end
if nargin < 9
  opts = [];
end

% validation phase
[objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts] = val.valptargin(objfun, x0, funvals0, true, lb, ub, lincon, nonlcon, multfun, opts);
opts.ValidateInput = false;
opts.OptSuppressOutput = true;

% initialization phase
[it0, it1, result, opts, stats] = pt.traceinit(objfun, lincon, nonlcon, x0, funvals0, opts);
[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(it1);
sizes = struct('variables', n, 'objectives', nobj, 'ineqlin', na, 'eqlin', naeq, 'ineqnonlin', nc, 'eqnonlin', nceq);

% stop condition (output function)
[initstop, ~, it1, ~, ~, stats] = pcoutfcn('INIT', [], it1, [], [], result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, []);
if initstop
  EXITFLAG = -1;
end

% first corrector
[~, it1, result, stats, stop] = pt.correct(it0, it1, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, @pcoutfcn);
if stop % stop condition (output function)
  EXITFLAG = 1;
  pt.tracefinit();
  return
end

if ~initstop
  if nobj == 2 && ~opts.PCForceCells
    % performs the continuation oriented to optimize the first function
    diropts = 'LeftUp';
    [it2, itlast, result, stats, EXITFLAG] = pt.trace2(it0, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, diropts, opts, stats, @pcoutfcn);
    if isempty(it2)
      it2 = it0;
    end
  
    % reverse the list
    result.ps(1 : stats.Count, :) = flipud(result.ps(1 : stats.Count, :));
    result.pf(1 : stats.Count, :) = flipud(result.pf(1 : stats.Count, :));
  
    if EXITFLAG == 1
      % performs the continuation oriented to optimize the second function
      diropts = 'RightDown';
      [~, itlast, result, stats, EXITFLAG] = pt.trace2(it2, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, diropts, opts, stats, @pcoutfcn);
    end
  else
    % performs the continuation for problems with more than two objectives.
    [itlast, result, stats, EXITFLAG] = pt.tracen(it0, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, @pcoutfcn);
  end
end

% finalization
[result, stats] = pt.tracefinit(result, stats);
[~, ~, ~, ~, ~, stats] = pcoutfcn('FINIT', [], itlast, [], [], result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG);

  function [stop, discard, it1, itp, itc, stats] = pcoutfcn(PHASEFLAG, it0, it1, itp, itc, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG)
    if ~isempty(opts.PCOutFcn)
      info = struct(...
        'PHASEFLAG', PHASEFLAG,...
        'it0', it0,... 
        'it1', it1,... 
        'itp', itp,...
        'itc', itc,...
        'objfun', objfun,... 
        'x0', x0,... 
        'funvals0', funvals0,... 
        'lb', lb,... 
        'ub', ub,... 
        'lincon', lincon,... 
        'nonlcon', nonlcon,... 
        'multfun', multfun,... 
        'opts', opts,... 
        'result', result,... 
        'stats', stats,... 
        'EXITFLAG', EXITFLAG,...
        'sizes', sizes);
      output = cell(1, nargout(opts.PCOutFcn));
      [output{:}] = feval(opts.PCOutFcn, info);
      l = length(output);
      if l > 0
        stop = output{1};
      else
        stop = false;
      end
      if l > 1
        discard = output{2};
      else
        discard = false;
      end
      if l > 2
        it1 = output{3};
      end
      if l > 3
        itp = output{4};
      end
      if l > 4
        itc = output{5};
      end
      if l > 5
        stats = output{6};
      end
    else
      discard = false;
      stop = false;
    end
  end
end

