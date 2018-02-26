% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [wJvx, wJvcount] = wjvtimes(wJx, v)
mx = size(wJx, 1);
mv = size(v, 1);
m = max(mx, mv);

wJvcount = m;
wJvx = sum(bsxfun(@times, wJx, v), 2); % (m x 1)
end

