function [G, m] = ugrid(n, m, lb, ub)
% Uniform n-grid of approximately m points uniformly distributed.
% If no bound is specified, each point component is between 0 and 1, i.e., 
% 0 <= x1 <= 1, ..., 0 <= xn <= 1.
% Returns G, a matrix of size (m x n).
% Returns the exact m required.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

mi = ceil(m^(1 / n)); % number of points for each dimension
m = mi^n;

if ~exist('lb', 'var') && ~exist('ub', 'var') ||...
   isempty(lb) && isempty(ub)
 
  xi = linspace(0, 1, mi);
  x = repmat({xi}, 1, n);
else
  if ~exist('lb', 'var') || ~exist('ub', 'var')
    error('utils:ugrid:LbOrUbMissing', 'Both lb and ub constraints must be specified.');
  end
  if isempty(lb) || isempty(ub)
    error('utils:ugrid:LbOrUbMissing', 'Both lb and ub constraints must be specified.');
  end
  [lb, ub] = val.valboxcon(lb, ub, n, true);
  
  x = cellfun(@(i) linspace(lb(i), ub(i), mi), num2cell(1 : n), 'UniformOutput', false);
end

G = cell(1, n);
[G{:}] = ndgrid(x{:});
G = cat(n + 1, G{:});
G = reshape(G, [m, n]);