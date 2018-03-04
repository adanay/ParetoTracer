function [value] = discM(x, a, A, B)
% A in {1, 2, ...}.
% B > 0.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

value = 1 - x(:, 1).^a .* cos(pi * A * x(:, 1).^B).^2;
end

