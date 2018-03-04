function [it1, stats, diropts, dirargout] = mindir_fmincon(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, opts, stats, diropts)
% Computes the search direction v for a minimization algorithm using the 
% fmincon algorithm to solve the subproblem below.
%  min   d
% (v,d)
%  s.t.  J*v + 1/2*v'*H*v - d <= 0,
%                         A*v <= -a + b
%                        Jc*v <= -c 
%                       Aeq*v  = -aeq + beq
%                      Jceq*v  = -ceq
%                  l - x <= v <= u - x
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.
 
% sizes
[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(it1);

% options for the direction subproblem
if nargin < 11 || isempty(diropts)
  diropts = optimset('Algorithm', 'sqp',...
    'GradObj', 'on', 'GradConstr', 'on',...
    'MaxIter', opts.OptDirMaxIts,...
    'Display', 'off');
end

% ensures active sets are computed
it1 = pt.active(it1, lb, ub, false, opts);

% evaluating the gradients
[it1, stats] = pt.jeval(it1, objfun, lb, ub, nonlcon, false, opts, stats);

% lower and upper bounds of the direction subproblem
% for (v, d)
dirlb = [lb(:) - it1.x(:); -Inf];
if na + naeq + nc + nceq == 0
  dirub = [ub(:) - it1.x(:); 0];
else
  dirub = [ub(:) - it1.x(:); Inf];
end

% linear inequalities 
%if strcmpi(opts.HessModif, 'off') || (strcmpi(opts.HessApprox, 'fd') && ~opts.FDForceHess) 
%  dirA = [lincon.A, zeros(na, 1); it1.Jcx, zeros(nc, 1); it1.Jx, zeros(nobj, 1)];
%  dirb = [-it1.ax(:); -it1.cx(:); zeros(nobj, 1)];
%else
  dirA = [lincon.A, zeros(na, 1); it1.Jcx, zeros(nc, 1)];
  dirb = [-it1.ax(:); -it1.cx(:)];
%end

% linear equalities
dirAeq = [lincon.Aeq, zeros(naeq, 1); it1.Jceqx, zeros(nceq, 1)]; 
dirbeq = [-it1.aeqx(:); -it1.ceqx(:)];  

% initial point for the direction subproblem
x0 = zeros(n + 1, 1); % non ascendent feasible direction
    
% using fmincon to solve the direction subproblem
[x, fval, EXITFLAG, output, lambda] = fmincon(...
  @dirfun, x0,... objective function and initial guess 
  dirA, dirb,... % linear inequalities 
  dirAeq, dirbeq,... % linear equalities 
  dirlb, dirub,... lower and upper bounds 
  @dirnonlcon,... % nonlinear constraints
  diropts);
    
if ~isempty(x)
  % direction and first order optimality measure
  x = x(:)';
  it1.v = x(1 : n);
  it1.v(it1.xLbActive & it1.v < 0) = 0;
  it1.v(it1.xUbActive & it1.v > 0) = 0;
  it1.d = x(n + 1);
  it1.FirstOrdOpt = abs(it1.d);
  
  % Lagrange multipliers
  it1.w.objectives = lambda.ineqnonlin(:)';
  
  lower = lambda.lower(:)';
  it1.w.lower = lower(1 : n);
  
  upper = lambda.upper(:)';
  it1.w.upper = upper(1 : n);
  
  ineqlin = lambda.ineqlin(:)';
  it1.w.ineqlin = ineqlin(1 : na);
  it1.w.ineqnonlin = ineqlin(na + 1 : na + nc);
  
  eqlin = lambda.eqlin(:)';
  it1.w.eqlin = eqlin(1 : naeq);
  it1.w.eqnonlin = eqlin(naeq + 1 : naeq + nceq);
end
    
stats.OptDirIts = stats.OptDirIts + output.iterations;
dirargout = {x, fval, EXITFLAG, output, lambda};
  
  function [fx, gx] = dirfun(x)
    % Objective function and gradient of the direction subproblem.
    fx = x(n + 1);
    if nargout > 1
      gx = [zeros(n, 1); 1];
    end
  end

  function [cx, ceqx, Jcx, Jceqx] = dirnonlcon(x)
    % Nonlinear inequality constraints of the direction subproblem.
    
    v = x(1 : n);
    d = x(n + 1);
    
    Jvx = it1.Jx * v;
    stats.JvCount = stats.JvCount + 1;
    
    [vHx, it1, stats] = pt.vheval(v(:)', it0, it1, objfun, lb, ub, multfun, opts, stats);
    vHx = shiftdim(vHx, 1);
    vHvx = vHx * v;
    
    % SD approach by default if vHvx is not positive. 
    % TODO: Can this be improved?
    normv = norm(v);
    if ~strcmp(opts.HessModif, 'off') && normv > 0
      for i = 1 : nobj
        vHxi = vHx(i, :);
        normvHxi = norm(vHxi);
        if normvHxi > 0 && vHvx(i) < 0
          vHvx(i) = (vHxi * vHxi') * normv / normvHxi;
          %vHvx(i) = v' * v;
        end
      end
    end
    
    cx = Jvx + 0.5 * vHvx - d * ones(nobj, 1);
    
    if nargout > 2
      Jcx = [it1.Jx + vHx, -ones(nobj, 1)]';
    end
    
    ceqx = [];
    Jceqx = [];
  end

%   function [W] = dirhessmult(~, lambda, v)
%     % Hessian multiply function of the direction subproblem.
%     
%     w = lambda.ineqnonlin;
%     v = v(1 : n);
%     
%     [Hwvx, it1, stats] = pt.hwveval(w(:)', v(:)', it0, it1, objfun, lb, ub, multfun, opts, stats);
%     W = [Hwvx'; 0];
%   end
 
%   function [W] = dirhess(~, lambda)
%     % Hessian multiply function of the direction subproblem.
%     
%     w = lambda.ineqnonlin;
%     
%     [Hwx, it1, stats] = pt.hweval(w, it0, it1, objfun, lb, ub, multfun, [], opts, stats);
%     W = [[Hwx, zeros(n, 1)]; zeros(1, n + 1)];
%   end
end
