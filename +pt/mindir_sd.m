% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [it1, stats, diropts, dirargout] = mindir_sd(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, opts, stats, diropts)
% Computes the search direction v for a minimization algorithm using the 
% quadprog algorithm to solve the subproblem below.
%  min   1/2*v'*v + d
% (v,d)
%  s.t.  
%          J*v <= d
%          A*v <= -a + b
%         Jc*v <= -c 
%        Aeq*v  = -aeq + beq
%       Jceq*v  = -ceq
%   l - x <= v <= u - x
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.
 
% sizes
[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(it1);

% options for the direction subproblem
if nargin < 11 || isempty(diropts)
  diropts = optimset('Algorithm', 'interior-point-convex',...
    'MaxIter', opts.OptDirMaxIts,...
    'Display', 'off');
end

% ensures active sets are computed
it1 = pt.active(it1, lb, ub, false, opts);

% evaluating the gradients
[it1, stats] = pt.jeval(it1, objfun, lb, ub, nonlcon, false, opts, stats);

% linear inequalities 
dirA = [it1.Jx, -ones(nobj, 1); lincon.A, zeros(na, 1); it1.Jcx, zeros(nc, 1)];
dirb = [zeros(nobj, 1); -it1.ax(:); -it1.cx(:)];

% adding active lower bound constraints
if it1.rLb > 0
  rA = zeros(it1.rLb, n); 
  rA(:, it1.xLbActive) = -1;
  dirA = [dirA; rA, zeros(it1.rLb, 1)];
  
  rb = it1.x(it1.xLbActive) - lb(it1.xLbActive);
  dirb = [dirb; rb'];
end

% adding active upper bound constraints
if it1.rUb > 0
  rA = zeros(it1.rUb, n); 
  rA(:, it1.xUbActive) = -1;
  dirA = [dirA; rA, zeros(it1.rUb, 1)];
  
  rb = -it1.x(it1.xUbActive) + ub(it1.xUbActive);
  dirb = [dirb; rb'];
end

% linear equalities
dirAeq = [lincon.Aeq, zeros(naeq, 1); it1.Jceqx, zeros(nceq, 1)]; 
dirbeq = [-it1.aeqx(:); -it1.ceqx(:)];  

% initial point for the direction subproblem
x0 = zeros(n + 1, 1); % non ascendent feasible direction
    
% using quadprog to solve the direction subproblem
[x, fval, EXITFLAG, output, lambda] = quadprog(...
  [speye(n), spalloc(n, 1, n); spalloc(1, n + 1, n + 1)], [spalloc(n, 1, n); 1],... subproblem function 
  ...[eye(n), zeros(n, 1); zeros(1, n + 1)], [zeros(n, 1); 1],...
  dirA, dirb,... inequality constraints
  dirAeq, dirbeq,... % linear equalities 
  [], [],...
  [],... starting point
  diropts); % options
    
if ~isempty(x)
  % direction and first order optimality measure
  x = x(:)';
  it1.v = x(1 : n);
  it1.v(it1.xLbActive & it1.v < 0) = 0;
  it1.v(it1.xUbActive & it1.v > 0) = 0;
  if na + naeq + nc + nceq == 0
    it1.v = it1.v / norm(it1.v);
  end
  it1.d = x(n + 1);
  it1.FirstOrdOpt = abs(it1.d);
  
  % Lagrange multipliers
  ineqlin = lambda.ineqlin(:)';
  it1.w.objectives = ineqlin(1 : nobj);
  it1.w.ineqlin = ineqlin(nobj + 1 : nobj + na);
  it1.w.ineqnonlin = ineqlin(nobj + na + 1 : nobj + na + nc);
  it1.w.lower = ineqlin(nobj + na + nc + 1 : nobj + na + nc + it1.rLb);
  it1.w.upper = ineqlin(nobj + na + nc + it1.rLb + 1 : nobj + na + nc + it1.rLb + it1.rUb);
  
  eqlin = lambda.eqlin(:)';
  it1.w.eqlin = eqlin(1 : naeq);
  it1.w.eqnonlin = eqlin(naeq + 1 : naeq + nceq);
end
    
stats.OptDirIts = stats.OptDirIts + output.iterations;
dirargout = {x, fval, EXITFLAG, output, lambda};
end
