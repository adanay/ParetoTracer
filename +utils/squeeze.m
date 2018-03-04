function [B] = squeeze(A, dim)
% Remove the specified singleton dimensions from A (not all as default).
% The algorithm just permutes the dimensions not wanted to the end so if 
% they are not singleton this approach does not work.

% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

ndimsA = ndims(A);
maxdim = max([ndimsA, dim]);

dim2keep = 1 : maxdim;
dim2keep(dim) = [];

B = permute(A, [dim2keep, dim]);
end

