% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [R] = choln(A, dim1, dim2, type)
% Computes the Cholesky factorization of all square matrices contained in 
% A(... : ... : ...).
%      dim1  dim2
% The necessary condition for the operation to be successful is that 
% size(A, dim1) = size(A, dim2) and that each A(... : ... : ...) is 
% positive definite.
% R is the Cholesky factorization matrix where 
% - size(R) = size(A) 
% - R(... : ... : ...) = chol(A(... : ... : ...))
%        dim1  dim2                dim1  dim2
% If type = 'lower', a lower triangular matrix is produced.

if nargin < 2
  dim1 = 1;
end
if nargin < 3
  dim2 = 2;
end
if nargin < 4
  type = 'upper';
end

if size(A, dim1) ~= size(A, dim2)
  error('utils:choln:ANotSquare', 'The size of the dimension %d of matrix A must coincide with the size of the dimension %d.', dim1, dim2);
end

dim = [dim1 dim2];
ndimsA = ndims(A);
maxdim = max([ndimsA, dim]);

R = cellfun(@(A) mycholn(A), num2cell(A, dim), 'UniformOutput', false);
R = cell2mat(R);

  function [R] = mycholn(A)
    dim2squeeze = 1 : maxdim; % squeeze all dimensions but dim
    dim2squeeze(dim) = [];
    A = permute(A, [dim, dim2squeeze]);
    
    R = chol(A, type);
    
    restore = 1 : maxdim; % restore permutation
    restore(dim2squeeze) = 3 : 3 + length(dim2squeeze) - 1;
    restore(dim) = [1, 2];
    
    R = permute(R, restore);
  end
end
