% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [value] = r_sum(y, i, w)
% w > 0.
% i represents the considered indices from y.

value = sum(bsxfun(@times, y(:, i), w), 2) / sum(w);
end

