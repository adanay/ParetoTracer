% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [A] = eyen(s, dim1, dim2)
% Returns a matrix A of the specified size such that
% A(... : ... : ...) = eye(n),
%      dim1  dim2
% where n = s(dim1) = s(dim2).

if nargin < 2
  dim1 = 1;
end
if nargin < 3
  dim2 = 2;
end

if s(dim1) ~= s(dim2)
  error('utils:eyen:ResultNotSquare', 'The size of the dimension %d must coincide with the size of the dimension %d.', dim1, dim2);
end

s = s(:)';
l = length(s);
n = s(dim1);

maxdim = max([dim1, dim2]);
if maxdim > l
  s = [s, ones(1, maxdim - l)];
end

z = ones(1, maxdim);
z(dim1) = n;
z(dim2) = n;
A = zeros(z);

i = cell(1, maxdim); 
i(:) = {1};
i(dim1) = {':'};
i(dim2) = {':'};
A(i{:}) = eye(n);

s(dim1) = 1;
s(dim2) = 1;
A = repmat(A, s);
end

