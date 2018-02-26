% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [JJx, JJcount] = jjtimes(Jx)
[m, nobj, n] = size(Jx);

JJcount = m;

%JJx = utils.mtimesn(permute(Jx, [1 3 4 2]), permute(Jx, [1 4 3 2]), 2, 3);

JJx(m, n, n, nobj) = 0;
for i = 1 : m
  for j = 1 : nobj
    Jxij = Jx(i, j, :); % (1 x 1 x n)
    Jxij = shiftdim(Jxij, 1);  % (1 x n)
    JJxij = Jxij' * Jxij; % (n x n)
    JJx(i, :, :, j) = JJxij;
  end
end
end

