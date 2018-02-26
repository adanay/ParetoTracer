% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Jvx, Jvcount] = jvtimes(Jx, v)
[mx, nobj, ~] = size(Jx);
[mv] = size(v, 1); 
m = max ([mx, mv]);

Jvcount = m;

%Jvx = utils.mtimesn(Jx, v, 2, 3);

Jvx(m, nobj) = 0;
for i = 1 : m
  ix = min(i, mx);
  iv = min(i, mv);
  Jxi = Jx(ix, :, :); % (1 x nobj x n)
  Jxi = shiftdim(Jxi, 1); % (nobj x n)
  Jvxi = Jxi * v(iv, :)'; % (nobj x 1)
  Jvx(i, :) = Jvxi';
end
end

