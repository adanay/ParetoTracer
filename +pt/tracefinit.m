function [result, stats] = tracefinit(result, stats)
% Fill the structures for the output of the PT algorithm. 

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

result.ps = result.ps(1 : stats.Count, :);
result.pf = result.pf(1 : stats.Count, :);

if stats.OptCount > 0
  stats.OptIts = stats.OptIts / stats.OptCount;
  stats.OptDirIts = stats.OptDirIts / stats.OptCount;
  stats.OptLsIts = stats.OptLsIts / stats.OptCount;
end
end

