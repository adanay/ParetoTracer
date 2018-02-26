% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [feasible] = isfeasible(x, lb, ub)
% Determines whether or not a point is feasible.
%
% x is the point to verify. For simplicity, x is assumed to be a row vector
% or a matrix of size (m x n) where m is the number of individuals and n is
% the number of variables.
% lb and ub are vectors that represent the box constraints.

lb = lb(:)';
ub = ub(:)';
feasible = all(bsxfun(@ge, x, lb), 2) & all(bsxfun(@le, x, ub), 2);
end



