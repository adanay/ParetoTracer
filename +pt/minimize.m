function [result, stats, EXITFLAG] = minimize(objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts)
% pt.minimize finds a constrained minimum of a function of several variables
% and several objectives.
%
% This version supports only equality constraints and box constrains.
% Support for general inequality constraints is still under construction.
% If some inequality constraint is passed, it is ignored.
%
% pt.minimize attempts to solve problems of the form:
%    min f(x)    subject to: A*x <= b, Aeq*x  = beq (linear constraints)
%     x                     c(x) <= 0, ceq(x) = 0   (nonlinear constraints)
%                             lb <= x <= ub         (box contraints)
% pt.minimize implements an algorithm where in each iteration a
% quadratically constrained linear subproblem is solved to compute a search
% direction that is coupled with a linear search strategy guided by the
% Armijo condition.
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
% x0: The initial guess. It must be a vector of n components.
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
% - x, fx, Jx, Hx: Respectively the solution, and the value of the
%   objective functions, Jacobian, and Hessians (if available).
% - ax and aeqx: Respectively the values of the linear inequality
%   and equality constraints (if available).
% - cx and ceqx: Respectively the values of the nonlinear inequality
%   and equality constraints (if available).
% - dcx and dceqx: Respectively the square norm of the nonlinear inequality
%   and equality constraints (if available).
% - Jcx, Jceqx: Respectively the values of the Jacobians of the nonlinear
%   inequality and equality constraints (if available).
% - w: Structure with the Lagrange multipliers.
%   + objectives: Objective functions.
%   + lower: Lower bounds.
%   + upper: Upper bounds.
%   + ineqlin: Linear inequalities.
%   + eqlin: Linear equalities.
%   + ineqnonlin: Nonlinear inequalities.
%   + eqnonlin: Nonlinear equalities.
% - v: Search direction (if computed) at the solution.
% - d: Measure of the objectives decrease.
% - t: Step length (if computed) at the solution.
% - FirstOrdOpt: Measure of the first-order optimality.
%
% EXITFLAG: Describes the exit condition.
%  0: Number of iterations exceeded opts.OptMaxIts.
%  1: First-order optimality measure was less than opts.FirstOrdOptTol.
%  2: Change in x too small.
% -1: Stopped by an output function.
% -2: No feasible point was found.
%
% stats: Statistics.
% - OptIts: Number of iterations taken. The iterations only increment if a 
%   new value of x is computed.
% - OptDirIts: Average iterations taken by the direction subproblem.
% - OptLsIts: Average backtrack iterations taken by the line search to
%   satisfy the Armijo condition.
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
[objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts] = val.valptargin(objfun, x0, funvals0, false, lb, ub, lincon, nonlcon, multfun, opts);
opts.ValidateInput = false; % do not validate again, e.g., by FD

% initialization phase
[it0, it1, stats] = pt.mininit(objfun, x0, funvals0, lincon, nonlcon, opts);
diropts = [];
[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(it1);
sizes = struct('variables', n, 'objectives', nobj, 'ineqlin', na, 'eqlin', naeq, 'ineqnonlin', nc, 'eqnonlin', nceq);

% stop condition (output function)
[initstop, it1, ~, stats] = optoutfcn('INIT', [], it1, [], [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats);
if initstop
  EXITFLAG = -1;
end

issd = pt.minissd(objfun, multfun, opts);

% iterations phase
% it0: previous iteration values
% it1: current iteration values
% it2: next iteration values
while true && ~initstop
  % stop condition (output function)
  [stop, it1, ~, stats] = optoutfcn('IT-INIT', it0, it1, [], [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats);
  if stop
    EXITFLAG = -1;
    [~, it1, ~, stats] = optoutfcn('IT-FINIT', it0, it1, [], [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG);
    break
  end
  
  % search direction
  if issd
    [it1, stats, diropts, dirargout] = pt.mindir_sd(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, opts, stats, diropts);
  else
    switch(lower(opts.OptDirSolver))
      case 'fmincon'
        [it1, stats, diropts, dirargout] = pt.mindir_fmincon(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, opts, stats, diropts);
    end
  end
  
  % stop condition (no feasible step found)
  if isempty(it1.v)
    EXITFLAG = -2; 
    [~, it1, ~, stats] = optoutfcn('IT-FINIT', it0, it1, [], [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG);
    break
  end
  
  % stop condition (first order optimality measure + constraint violation measure) 
  if it1.FirstOrdOpt <= opts.FirstOrdOptTol && it1.dceqx <= opts.ConstViolTol;
    EXITFLAG = 1;
    [~, it1, ~, stats] = optoutfcn('IT-FINIT', it0, it1, [], [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG);
    break
  end
  
  % stop condition (output function)
  [stop, it1, ~, stats] = optoutfcn('DIR', it0, it1, [], [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, [], dirargout);
  if stop
    EXITFLAG = -1;
    [~, it1, ~, stats] = optoutfcn('IT-FINIT', it0, it1, [], [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG);
    break
  end
  
  % step length (line search)
  [it1, it2, stats, lsargout] = pt.minls(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, opts, stats);
  
  % stop condition (output function)
  [stop, it1, it2, stats] = optoutfcn('STEP', it0, it1, it2, [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, [], [], lsargout);
  if stop
    EXITFLAG = -1;
    [~, it1, ~, stats] = optoutfcn('IT-FINIT', it0, it1, it2, [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG);
    break
  end
  
  % stop condition (change in x too small)
  if it1.t == 0; % actually no step length was found
    EXITFLAG = 2;
    [~, it1, ~, stats] = optoutfcn('IT-FINIT', it0, it1, it2, [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG);
    break
  end
  
  % stop condition (change in x too small)
  if norm(it2.x - it1.x) <= opts.OptMinStepVar;
    EXITFLAG = 2;
    [~, it1, ~, stats] = optoutfcn('IT-FINIT', it0, it1, it2, [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG);
    break
  end
  
  [stop, it1, it2, stats] = optoutfcn('IT-FINIT', it0, it1, it2, [], objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats);
  if stop
    EXITFLAG = -1;
    break
  end
  
  % the iterations only increment if a new value of x is computed,
  % i.e., the algorithm did not stop before
  stats.OptIts = stats.OptIts + 1;
  
  % stop condition (max number of iterations reached)
  if stats.OptIts == opts.OptMaxIts
    EXITFLAG = 0;
    break
  end
  
  % reset
  it0 = it1;
  it1 = it2;
end

% finalization phase
[result, stats] = pt.minfinit(it1, stats);
[~, ~, ~, stats] = optoutfcn('FINIT', [], [], [], result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG);

  function [stop, it1, it2, stats] = optoutfcn(PHASEFLAG, it0, it1, it2, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, EXITFLAG, dirargout, lsargout)
    
    if ~exist('EXITFLAG', 'var')
      EXITFLAG = [];
    end

    if ~exist('dirargout', 'var') || isempty(dirargout)
      dirargout = [];
      DIREXITFLAG = [];
      OptDirIts = [];
    else
      switch(lower(opts.OptDirSolver))
        case 'fmincon'
          dirargout = struct(...
            'x', dirargout{1},...
            'fval', dirargout{2},...
            'exitflag', dirargout{3},...
            'output', dirargout{4},...
            'lambda', dirargout{3});
          DIREXITFLAG = dirargout.exitflag;
          OptDirIts = dirargout.output.iterations;
      end
    end
    
    if ~exist('lsargout', 'var') || isempty(lsargout)
      OptLsIts = [];
      LSEXITFLAG = [];
    else
      OptLsIts = lsargout{1};
      LSEXITFLAG = lsargout{2};
    end
    
    % info is a struct that gathers all current variables being used by the 
    % algorithm
    if ~isempty(opts.OptOutFcn)
      info = struct(...
        'PHASEFLAG', PHASEFLAG,...
        'it0', it0,...
        'it1', it1,...
        'it2', it2,...
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
        'sizes', sizes,...
        'dirargout', dirargout,...
        'DIREXITFLAG', DIREXITFLAG,...
        'OptDirIts', OptDirIts,...
        'LSEXITFLAG', LSEXITFLAG,...
        'OptLsIts', OptLsIts);

      % calling the output function
      output = cell(1, nargout(opts.OptOutFcn));
      [output{:}] = feval(opts.OptOutFcn, info);
      
      % the output function can replace the following variables
      l = length(output);
      if l > 0
        stop = output{1};
      else
        stop = false;
      end
      if l > 1
        it1 = output{2};
      end
      if l > 2
        it2 = output{3};
      end
      if l > 3
        stats = output{4};
      end
    else % no output function provided
      stop = false;
    end
  end
end


