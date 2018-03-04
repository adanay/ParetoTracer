%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.
function [Hwvx, Hwvcount] = hwvtimes(Hwx, v)
[mx, n, ~] = size(Hwx); 
[mv] = size(v, 1); 
m = max ([mx, mv]);

Hwvcount = m;

% Hwvx = utils.mtimesn(Hwx, v, 2, 3);
 
Hwvx(m, n) = 0;
for i = 1 : m
  ix = min(i, mx);
  iv = min(i, mv);
  Hwxi = Hwx(ix, :, :); % (1 x n x n)
  Hwxi = shiftdim(Hwxi, 1); % (n x n)
  Hwvxi = v(iv, :) * Hwxi; % (1 x n)
  Hwvx(i, :) = Hwvxi;
end
end