function [D] = detn(A, dim1, dim2)
% Computes the determinants of all square matrices contained in A(... dim1 ... dim2 ...).
% The necessary condition for the operation to be successful is that 
% size(A, dim1) = size(A, dim2).
% D is the determinat matrix where 
% - ndims(D) = ndims(A) 
% - size(D, dim1) = 1
% - size(D, dim2) = 1
% - size(D, dim?) = size(A, dim?)
% - D(... dim1 ... dim2 ...) = det(A(... dim1 ... dim2 ...))

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if nargin < 2
  dim1 = 1;
end
if nargin < 3
  dim2 = 2;
end

if size(A, dim1) ~= size(A, dim2)
  error('utils:detn:ANotSquare', 'The size of the dimension %d of matrix A must coincide with the size of the dimension %d.', dim1, dim2);
end

dim = [dim1 dim2];
D = cellfun(@(A) det(utils.slicen1(A, dim)), num2cell(A, dim));
end

