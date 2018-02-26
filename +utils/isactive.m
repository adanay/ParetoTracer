% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [a, la, ua] = isactive(x, lb, ub, tol)
% Active set using logical indices.
% lb and ub are the box constraints.
% x is the point to verify. For simplicity, x is assumed to be a row vector
% or a matrix of size (m x n) where m is the number of individuals and n is
% the number of variables.
% tol > 0 is a tolerance. 

[m, n] = size(x); % m is the number of individuals and n is the number of variables.

% active lower bound constraints
if ~isempty(lb)
  la = bsxfun(@ge, bsxfun(@plus, -x, lb), -tol);
else
  la = false(m, n);
end

% active upper bound constraints
if ~isempty(ub)
  ua = bsxfun(@ge, bsxfun(@minus, x, ub), -tol);
else
  ua = false(m, n);
end

% logical indices corresponding to active constraints
a = la | ua;
end
