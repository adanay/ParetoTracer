% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Hwx, Hwcount] = hwtimes(Hx, w)
mx = size(Hx, 1);
mw = size(w, 1); 
m = max (mx, mw);

Hwcount = m;

Hwx = sum(bsxfun(@times, Hx, permute(w, [1 3 4 2])), 4); % (m x n x n)
end

