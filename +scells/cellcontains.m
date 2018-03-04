function [contained] = cellcontains(c, r, x, lb, ub)
% Determines whether or not a point belongs to a cell.
%
% c is the center of the cell.
% r is the radius of the cell.
% x is the point to verify. For simplicity, x is assumed to be a row vector
% or a matrix of size (m x n) where m is the number of individuals and n is
% the number of variables.
% lb and ub are vectors that represent the box constraints.
%
% Returns true if the point belongs to the cell or false otherwise.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if isempty(c) || isempty(r)
  contained = false(m, 1);
  return
end

c = c(:)';
r = r(:)';
lb = lb(:)';
ub = ub(:)';

contained = all(bsxfun(@ge, x, c - r), 2) &... % points on the lower bound belong to the cell 
            all(bsxfun(@lt, x, c + r) | bsxfun(@eq, x, c + r) & bsxfun(@eq, x, ub), 2); 
            % points on the upper bound do not belong to the cell unless this is the upper bound of the space
end

