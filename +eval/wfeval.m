function [wfx, wfcount, fx, fcount, fundef] = wfeval(f, x, fx, w, iswarning, opts)
% Vectorized objective multiply function evaluation.
% wf = w' * f.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% x and w must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% function value is undefined. The function value will be checked only if 
% opts.FunValCheck is true.
% opts are the optimization options.
%
% The result is always of size (m x 1).
% The result of wJ(x) is assumed to be (1 x n) if m = 1.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

fcount = 0;
fundef = false;

mx = size(x, 1); % m is the number of individuals.
mw = size(w, 1); 
m = max(mx, mw);

if isempty(fx)
  [fx, fcount, fundef] = eval.feval(f, x, iswarning, opts);
end
wfx = sum(bsxfun(@times, w, fx), 2); % (m x 1)
wfcount = m;
end