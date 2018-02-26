% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [b] = cellneighbors(c, r, lb, ub)
% Determines the neighboring cells in coordinate alignment.
%
% c and r are the center and radius of the specified cell.
% lb and ub are vectors that represent the box constraints.
%
% Returns the neighboring cells that is a list of the respective centers.

c = c(:)';
r = r(:)';
lb = lb(:)';
ub = ub(:)';

n = length(c); % number of variables
m = 0; % number of neighbors
b = zeros(2 * n, n); % set of neighbors

for i = 1 : n
  % centers of the neighbor cells in the i-th dimension
  c1 = c;
  c2 = c;
    
  c1(i) = c1(i) - 2 * r(i);
  c2(i) = c2(i) + 2 * r(i);
    
  % feasibility of the left neighbor
  if c1(i) > lb(i)
    m = m + 1;
    b(m, :) = c1;
  end
    
  % feasibility of the right neighbor
  if c2(i) < ub(i)
    m = m + 1;
    b(m, :) = c2;
  end
end

b = b(1 : m, :);
end

