% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [R1, p] = cholupdaten(R, X, dim1, dim2, o)
% Updates the Cholesky factorization of all square matrices contained in
% R(... : ... : ...) with X(... : ... : ...).
%      dim1  dim2              dim1  dim2
% The necessary condition for the operation to be successful is that
% - the sizes of both R and X are identical except for size(X, dim2) = 1.
% - size(R, dim1) = size(R, dim2).
% o is either '+' or '-'.
%
% R1 is the updated Cholesky factor where
% - size(R1) = size(R)
% - R1(... : ... : ...) = cholupdate(R(... : ... : ...), X(... : ... 1 ...), o)
%         dim1  dim2                      dim1  dim2          dim1  dim2

if nargin < 3
  dim1 = 1;
end
if nargin < 4
  dim2 = 2;
end
if nargin < 5
  o = '+';
end

if size(R, dim1) ~= size(R, dim2)
  error('utils:cholupdaten:RNotSquare', 'The size of the dimension %d of matrix R must coincide with the size of the dimension %d.', dim1, dim2);
end

dim = [dim1 dim2];
ndimsR = ndims(R);
maxdim = max([ndimsR, dim]);

[R1, p] = cellfun(@(R, X) mycholupdaten(R, X), num2cell(R, dim), num2cell(X, dim), 'UniformOutput', false);
R1 = cell2mat(R1);
p = cell2mat(p);

  function [R1, p] = mycholupdaten(R, X)
    dim2squeeze = 1 : maxdim; % squeeze all dimensions but dim
    dim2squeeze(dim) = [];
    R = permute(R, [dim, dim2squeeze]);
    X = permute(X, [dim, dim2squeeze]);
    
    if o == '+'
      R1 = cholupdate(R, X, o);
      p = [];
    else
      [R1, p] = cholupdate(R, X, o);
    end
    
    restore = 1 : maxdim; % restore permutation
    restore(dim2squeeze) = 3 : 3 + length(dim2squeeze) - 1;
    restore(dim) = [1, 2];
    
    R1 = permute(R1, restore);
   
    if o ~= '+'
      p = permute(p, restore);
    end
  end
end

