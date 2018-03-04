%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.
function [Hwvx, Hwvcount] = vhwtimes(vHx, w)
mx = size(vHx, 1); 
mw = size(w, 1);
m = max (mx, mw);

Hwvcount = m;

Hwvx = sum(bsxfun(@times, vHx, w), 2); % (m x 1 x n)
Hwvx = permute(Hwvx, [1 3 2]); % (m x n)   
end