% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = s_decept(y, A, B, C)
% A in (0, 1).
% 0 < B << 1.
% 0 < C << 1.
% A - B > 0.
% A + B < 1.

t = abs(y - A) - B;
t1 = floor(y - A + B) * (1 - C + (A - B) / B) / (A - B);
t2 = floor(A + B - y) * (1 - C + (1 - A - B) / B) / (1 - A - B);

value = 1 + t .* (t1 + t2 + 1 / B);
end

