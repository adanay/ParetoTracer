function [L, modif] = modcholn(A, dim1, dim2)
% Computes the modified Cholesky factorization of all square matrices 
% contained in A(... : ... : ...).
%                   dim1  dim2
% The necessary condition for the operation to be successful is that 
% size(A, dim1) = size(A, dim2) and that each A(... : ... : ...) is
% positive definite.
% L is the Cholesky factorization matrix where 
% - size(L) = size(A) 
% - L(... : ... : ...) = modchol(A(... : ... : ...))
%        dim1  dim2                   dim1  dim2
% - modif(... : ... : ...) = true/false (whether the matrix was modified). 
%           dim1  dim2    

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
  error('utils:modcholn:ANotSquare', 'The size of the dimension %d of matrix A must coincide with the size of the dimension %d.', dim1, dim2);
end

dim = [dim1 dim2];
ndimsA = ndims(A);
maxdim = max([ndimsA, dim]);

[L, modif] = cellfun(@(A) mymodcholn(A), num2cell(A, dim), 'UniformOutput', false);
L = cell2mat(L);
modif = cell2mat(modif);

  function [L, modif] = mymodcholn(A)
    dim2squeeze = 1 : maxdim; % squeeze all dimensions but dim
    dim2squeeze(dim) = [];
    
    A = permute(A, [dim, dim2squeeze]);
    
    [L, modif] = utils.modchol(A);  
    
    restore = 1 : maxdim; % restore permutation
    restore(dim2squeeze) = 3 : 3 + length(dim2squeeze) - 1;
    restore(dim) = [1, 2];
    
    L = permute(L, restore);
    modif = permute(modif, restore);
  end
end
