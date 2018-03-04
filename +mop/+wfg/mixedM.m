%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = mixedM(x, a, A)
% a > 0.
% A = {1, 2, 3,...}.
value = (1 - x(:, 1) - cos(2 * A * pi * x(:, 1) + pi / 2) / (2 * A * pi)).^a;
end

