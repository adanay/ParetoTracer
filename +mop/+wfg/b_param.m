% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = b_param(y, i, j, w, A, B, C)
% w > 0.
% A in (0, 1).
% 0 < B < C.
% i represents the considered index from y.
% j represents the considered indices from y to be used by r_sum.
% i is not in the set of indices j.

u = mop.wfg.r_sum(y, j, w);
v = B + (C - B) * (A - (1 - 2 * u) .* abs(floor(0.5 - u) + A));

value = y(:, i).^v;
end

