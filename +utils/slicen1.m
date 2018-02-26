% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [S] = slicen1(A, dim)
% Takes an n-dimensional slice from A such that 
% S = A(1, 1, ..., :, ..., :, ..., :, 1,...1) 
%                dim1 ...dim2 ...dimn
% dim = [dim1, dim2, ... dimn].
% S is the slice matrix where ndims(S) = length(dim). 

ndimsA = ndims(A);
maxdim = max([ndimsA, dim]);

ind = cell(1, maxdim); 
ind(:) = {1};
ind(dim) = {':'};

S = A(ind{:});

dim2squeeze = 1 : maxdim;
dim2squeeze(dim) = [];
S = utils.squeeze(S, dim2squeeze);
end
