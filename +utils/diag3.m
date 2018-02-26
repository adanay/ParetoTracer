% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [D] = diag3(A)
% 3-D diagonal.
% If A is a matrix of size (m x n), D will be a 3-D matrix of size 
% (m x n x n) where D(i, :, :) is a diagonal matrix formed with the i-th 
% row of A.
% If A is a 3-D matrix of size(m x n x n), D will be a 2-D matrix of size 
% (m x n) where the i-th row of D is the diagonal of A(i, :, :).
if ismatrix(A)
  [m, n] = size(A);
  D = num2cell(A, 1);
  D = blkdiag(D{:});
  D = reshape(D, m, n, n);
else
  A = A(:, :, :);
  [m, n, p] = size(A);
  D = A(repmat(shiftdim(logical(eye(n, p)), -1), [m 1 1]));
  D = reshape(D, m, n);
end
end

