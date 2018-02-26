% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [dominant] = isdominant(y1, y2)
% Defines the Pareto optimality regarding minimization.
% y1 and y2 are vectors.
% For simplicity, they must be row vectors of size (m x n) where m is the
% number of individuals and n is the number of variables.
% Returns true if y1 dominates y2 and false otherwise.

dominant = all(bsxfun(@le, y1, y2), 2) && any(bsxfun(@lt, y1, y2), 2);
end

