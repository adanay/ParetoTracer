% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = b_flat(y, A, B, C)
% A, B, C in [0, 1].
% B < C.
% B = 0 => A = 0 and C ~= 1.
% C = 1 => A = 1 and B ~= 0.

u =       A * min(0, floor(y - B)) .* (B - y) / B;
v = (1 - A) * min(0, floor(C - y)) .* (y - C) / (1 - C);
value = A + u - v;
value = roundn(value, -6);
end

