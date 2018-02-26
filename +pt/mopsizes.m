% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [n, nobj, na, naeq, nc, nceq] = mopsizes(it1)
% Returns the sizes of a MOP.
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.

n = length(it1.x);
nobj = length(it1.fx);
na = length(it1.ax);
naeq = length(it1.aeqx);
nc = length(it1.cx);
nceq = length(it1.ceqx);
end

