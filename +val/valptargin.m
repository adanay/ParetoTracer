function [objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts] = valptargin(objfun, x0, funvals0, checkW, lb, ub, lincon, nonlcon, multfun, opts)
% Validates the arguments of a Pareto Tracer algorithm.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

doval = ~(nargin >= 7 && isfield(opts, 'ValidateInput') && isa(opts.ValidateInput, 'logical') && ~any(opts.ValidateInput));

% objective function validation
objfun = val.valobjfun(objfun, doval);

% starting point validation
x0 = x0(:)';
n = length(x0);
if isempty(x0)
  error('val:valptargin:EmptyX0', 'An initial guess x0 must be specified.')
end
val.checkval(x0, 'x0', false);

% options validations
%jacprovided = isstruct(info.objfun) && ~isempty(info.objfun.J);
%hessprovided = isstruct(info.objfun) && ~isempty(info.objfun.H);
opts = val.valptopts(opts, n, doval);

% box constraints validation
[lb, ub] = val.valboxcon(lb, ub, n, doval);
if doval && any(x0 < lb) || any(x0 > ub)
  error('val:valptargin:InfeasX0', 'The initial point x0 is infeasible.')
end

% linear constraints validation
[lincon, na, naeq] = val.vallincon(lincon, n, doval);

% starting function values validation
funvals0 = val.valfunvals(funvals0, checkW, n, na, naeq, opts.FunValCheck, doval);

% nonlinear constraints validation
nonlcon = val.valnonlcon(nonlcon, doval);

% multiply functions validation
multfun = val.valhessmultfun(multfun, doval);
end