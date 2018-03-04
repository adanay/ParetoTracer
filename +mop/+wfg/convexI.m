function [value] = convexI(x)
% x is supposed to have 2 or more components.

% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

value = fliplr(cumprod(1 - cos(x(:, 1 : end - 1)  * pi / 2), 2)) .*...
                      (1 - sin(x(:, end : -1 : 2) * pi / 2));
end

