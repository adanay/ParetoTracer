function [value] = r_sum(y, i, w)
% w > 0.
% i represents the considered indices from y.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

value = sum(bsxfun(@times, y(:, i), w), 2) / sum(w);
end

