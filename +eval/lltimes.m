% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hx, LLcount] = lltimes(Lx)
[m, n, ~, nobj] = size(Lx);

LLcount = m;

%Hx = utils.mtimesn(Lx, permute(Lx, [1 3 2 4]), 2, 3);

Hx(m, n, n, nobj) = 0;
for i = 1 : m
  for j = 1 : nobj
    Lxij = Lx(i, :, :, j); % (1 x n x n)
    Lxij = shiftdim(Lxij, 1);  % (n x n)
    Hxij = Lxij * Lxij'; % (n x n)
    Hx(i, :, :, j) = Hxij;
  end
end
end