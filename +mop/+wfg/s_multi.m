% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = s_multi(y, A, B, C)
% A in {1, 2,...}.
% B >= 0.
% (4A + 2)pi >= 4B.
% C in (0, 1).

r = floor(C - y) + C;
t = abs(y - C) ./ (2 * r);
s = (4 * A + 2) * pi * (0.5 - t);
z = 4 * B * t.^2;

value = (1 + cos(s) + z) / (B + 2);
end

