% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = concaveI(x)
% x is supposed to have 2 or more components.

value = fliplr(cumprod(sin(x(:, 1 : end - 1)  * pi / 2), 2)) .*...
                       cos(x(:, end : -1 : 2) * pi / 2);
end

