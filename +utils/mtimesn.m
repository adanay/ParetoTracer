% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function C = mtimesn(A, B, dim1, dim2)
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

if nargin < 3
  dim1 = 1;
end
if nargin < 4
  dim2 = 2;
end

if(ismatrix(A) && ismatrix(B) && dim1 == 1 && dim2 == 2)
  C = A * B;
  return
end

if size(A, dim2) ~= size(B, dim1)
  error('vec:mtimesn:ABNotConsistent', 'The size of the dimension %d of matrix A must coincide with the size of the dimension %d of matrix B.', dim2, dim1);
end

% TODO: This is very slow. Find alternatives. 
C = utils.zzz_mmat(A, B, dim1, dim2);
end

