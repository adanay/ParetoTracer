%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.
function [vHx, vHcount] = vhtimes(Hx, v)
[mx, n, ~, nobj] = size(Hx); 
mv = size(v, 1); 
m = max(mx, mv);

vHcount = m;

% vHx = utils.mtimesn(permute(v, [1 3 2]), Hx, 2, 3);
% vHx = permute(vHx, [1 4 3 2]);

vHx(m, nobj, n) = 0;
for i = 1 : m
  ix = min(i, mx);
  iv = min(i, mv);
  for j = 1 : nobj
    Hxij = Hx(ix, :, :, j); % (1 x n x n)
    Hxij = shiftdim(Hxij, 1); % (n x n)
    vHxij = v(iv, :) * Hxij; % (1 x n)
    vHx(i, j, :) = vHxij;
  end
end
end
