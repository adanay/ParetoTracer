% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function C = zzz_mmat(A, B, dim1, dim2)
% Simple matrix multiplication of multidimensional arrays.
%
% Input:
% A, B        Multidimensional input arrays.
% dim1, dim2  The two dimensions of multiplication.
%
% The two matrices are selected by:
% C = A(... dim1 ... dim2 ...) * B(... dim1 ... dim2 ...)
% The necessary condition for the multiplication to be sucessfully 
% performed is that size(A, dim2) = size(B, dim1).
% Singleton dimensions in both A and B matrices are supported.
% The default values for dim1 and dim2 is 1 and 2 respectively.
%
% Examples:
% For 2D matrices mmat is identical to the Matlab built-in multiplication:
% A = [1 2 3 4];
% B = [1;2;3;4];
% C = mmat(A,B)
%
% C will be 30.
%
% For multidimensional arrays:
% A = repmat([1 2 3 4],[1 1 5]);
% B = [1 2 3 4]';
% C = mmat(A,B)
% C will be an array with dimensions of 1x1x5 and every element is 30.

dim = [dim1 dim2];
ndimsA = ndims(A);
ndimsB = ndims(B);
maxdim = max([ndimsA, ndimsB, dim]);

sizeA = [size(A), ones(1, maxdim - ndimsA)]; 
sizeA = sizeA(dim);
sizeB = [size(B), ones(1, maxdim - ndimsB)]; 
sizeB = sizeB(dim);

% form A matrix
A = repmat(A, [ones(1, maxdim) sizeB(2)]);

% form B matrix
ind = 1 : maxdim + 1; 
ind([dim end]) = ind([end dim]);
rep = ones(1, maxdim + 1); 
rep(dim(1)) = sizeA(1);
B = repmat(permute(B, ind), rep);

% multiply with expanding along singleton dimensions
C = sum(bsxfun(@times, A, B), dim2);

% form C matrix
ind = 1 : maxdim + 1; 
ind([dim end]) = ind([dim(1) end dim(2)]);
C = permute(C, ind);
end