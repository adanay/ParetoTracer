function [x] = project(x, lb, ub)
% Project the value of x onto the box. For simplicity, x is assumed to be a 
% row vector or a matrix of size (m x n) where m is the number of 
% individuals and n is the number of variables.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

m = size(x, 1);
if ~isempty(lb)
  lb = repmat(lb, m, 1);
  i = x < lb;
  x(i) = lb(i);
end
if ~isempty(ub)
  ub = repmat(ub, m, 1);
  i = x > ub;
  x(i) = ub(i);
end
end
