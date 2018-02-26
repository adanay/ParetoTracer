% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [B] = squeeze(A, dim)
% Remove the specified singleton dimensions from A (not all as default).
% The algorithm just permutes the dimensions not wanted to the end so if 
% they are not singleton this approach does not work.

ndimsA = ndims(A);
maxdim = max([ndimsA, dim]);

dim2keep = 1 : maxdim;
dim2keep(dim) = [];

B = permute(A, [dim2keep, dim]);
end

