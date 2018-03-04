function [it1, result, opts, stats] = traceresult(it1, opts, stats)
% Initializes the Pareto Tracer result structure. 
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% sizes
n = length(it1.x);
nobj = length(it1.fx);
N = opts.InitSetSize;

% Pareto set and front
ps(N, n) = 0;
pf(N, nobj) = 0;

result = struct(...
  'ps', ps,...
  'pf', pf,...
  'tree', [],...
  'depth', [],...
  'cellradius', []);

% this validation cannot take place until the function is first evaluated
doval = ~(isfield(opts, 'ValidateInput') && isa(opts.ValidateInput, 'logical') && ~any(opts.ValidateInput));
[lbobj, ubobj] = val.valboxcon(opts.LbObj, opts.UbObj, nobj, doval);
opts.LbObj = lbobj;
opts.UbObj = ubobj;

if nobj > 2 || opts.PCForceCells
  result.depth = nobj * ceil(log2(max(ubobj - lbobj) / opts.PCStepObj));
  result.tree = scells.treenode;
  result.cellradius = scells.radius(lbobj, ubobj, result.depth);
end
end

