% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = linearI(x)
% x is supposed to have 2 or more components.
% Supposing M = 5,
% value = x1 * x2 * x3 * (1 - x4)
%         x1 * x2 * (1 - x3)
%         x1 * (1 - x2)

value = fliplr(cumprod(x(:, 1 : end - 1), 2)) .*...
                  (1 - x(:, end : -1 : 2));
end

