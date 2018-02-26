% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [wJvx, wJvcount] = jvwtimes(Jvx, w)
mx = size(Jvx, 1);
mw = size(w, 1);
m = max (mx, mw);

wJvcount = m;
wJvx = sum(bsxfun(@times, Jvx, w), 2); % (m x 1)
end

