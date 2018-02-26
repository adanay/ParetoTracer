% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [N] = normn(A, dim1, dim2, p)
% Computes the p-norm of all matrices contained in A(... dim1 ... dim2 ...).
% N is the p-norm matrix where 
% - ndims(N) = ndims(A) 
% - size(N, dim1) = 1
% - size(N, dim2) = 1
% - size(N, dim?) = size(A, dim?)
% - N(... dim1 ... dim2 ...) = norm(A(... dim1 ... dim2 ...), p)

if nargin < 2
  dim1 = 1;
end
if nargin < 3
  dim2 = 2;
end
if nargin < 4
  p = 2;
end

dim = [dim1 dim2];
N = cellfun(@(A) norm(utils.slicen1(A, dim), p), num2cell(A, dim));
end
