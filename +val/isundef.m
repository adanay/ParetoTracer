function [undef] = isundef(M)
% Determines whether a matrix has some undefined entry.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

undef = ~isreal(M) || any(isinf(M(:))) || any(isnan(M(:)));
end

