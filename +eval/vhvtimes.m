% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [vHvx, vHvcount] = vhvtimes(vHx, v)
[mx, nobj, ~] = size(vHx); 
mv = size(v, 1);
m = max (mx, mv);

vHvcount = m;

% vHvx = utils.mtimesn(vHx, v, 2, 3);

vHvx(m, nobj) = 0;
for i = 1 : m
  ix = min(i, mx);
  iv = min(i, mv);
  for j = 1 : nobj
    vHxij = vHx(ix, j, :); % (1 x 1 x n)
    vHxij = shiftdim(vHxij, 2); % (n x 1)
    vHvxij = v(iv, :) * vHxij; % (1 x n)
    vHvx(i, j, :) = vHvxij;
  end
end   
end