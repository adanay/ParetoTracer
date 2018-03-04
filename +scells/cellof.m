function [c, r] = cellof(x, lb, ub, depth)
% Determines the cell covering the entered point.
%
% x is the point whose covering cell needs to be determined.
% lb and ub are vectors that represent the box constraints.
% depth is the number of subdivision iterations.
%
% Returns the center and the radius of the cell.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

x = x(:)';
lb = lb(:)';
ub = ub(:)';
n = length(x); % number of variables

% feasibility of the point
if ~utils.isfeasible(x, lb, ub)
  c = [];
  r = [];
  return;
end

r = scells.radius(lb, ub, depth); % radius

% lower bound
l = lb + 2 * r .* max(0, floor((x - lb) ./ (2 * r))); 
for i = 1 : n
  if l(i) >= ub(i)
    l(i) = l(i) - 2 * r(i);
  end
end

c = l + r; % center
end
