% 

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.
function [wJx, wJcount] = wjtimes(Jx, w)
mx = size(Jx, 1);
mw = size(w, 1);
m = max(mx, mw);

wJcount = m;

wJx = sum(bsxfun(@times, w, Jx), 2); % (m x 1 x n)
wJx = permute(wJx, [1, 3, 2]); % (m x n)
end

