% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = s_linear(y, A)
% A in (0, 1).

value = abs(y - A) ./ abs(floor(A - y) + A);
end

