function [value] = s_linear(y, A)
% A in (0, 1).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

value = abs(y - A) ./ abs(floor(A - y) + A);
end

