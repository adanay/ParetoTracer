% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [itc, result, stats] = tracesave(itc, opts, result, stats)
% sizes
n = length(itc.x);
nobj = length(itc.fx);

% resizing sets
if stats.Count == size(result.ps, 1)
  N = ceil(stats.Count * opts.SetSizeGrowFact) + 1;
  result.ps(N, n) = 0;
  result.pf(N, nobj) = 0;
end

stats.Count = stats.Count + 1;

% saving the corrector
result.ps(stats.Count, :) = itc.x;
result.pf(stats.Count, :) = itc.fx;
itc.index = stats.Count;

% update the tree
if ~isempty(result.tree)
  [c, ~, inserted, index] = scells.recover(result.tree, result.depth, itc.fx, opts.LbObj, opts.UbObj, stats.Count);
  itc.cell = c;
  itc.cellinserted = inserted;
  itc.index = index;
end
end
