% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [r] = radius(lb, ub, depth)
% Determines the radius of cells at the current subdivision iteration.
%
% lb and ub are vectors that represent the box constraints.
% depth is the number of subdivision iterations.
%
% Returns the radius of cells. 

n = length(lb); % number of variables
j = mod(depth - 1, n) + 1; % current subdivision coordinate

r = zeros(1, n);
for i = 1 : n
  if i <= j
    r(i) = (ub(i) - lb(i)) / 2^(floor((depth - 1) / n) + 2);
  else
    r(i) = (ub(i) - lb(i)) / 2^(floor((depth - 1) / n) + 1);
  end
end
end

