%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.
function [vHx, vHcount] = vlltimes(Lx, v)
[mx, n, ~, nobj] = size(Lx);
mv = size(v, 1); 
m = max(mx, mv); 

vHcount = m;

vHx(m, nobj, n) = 0;
for i = 1 : m
  ix = min(i, mx);
  iv = min(i, mv);
  for j = 1 : nobj
    Lxij = Lx(ix, :, :, j); % (1 x n x n)
    Lxij = shiftdim(Lxij, 1);  % (n x n)
    vHxij = v(iv, :) * Lxij; % (1 x n)
    vHxij = vHxij * Lxij'; % (1 x n)
    vHx(i, j, :) = vHxij;
  end
end
end