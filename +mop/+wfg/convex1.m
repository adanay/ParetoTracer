%

% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = convex1(x)
value = prod(1 - cos(pi * x / 2), 2);
end

