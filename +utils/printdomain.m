function [] = printdomain(lb, ub, fid)
% Prints the string representation of a domain set.
% lb and ub are vectors that represent the box constraints.
% fid is an optional parameter that represents a file id.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

n = length(lb); % number of variables

for i = 1 : n
  if exist('fid', 'var')
    fprintf(fid, '[%.2f,%.2f]', lb(i), ub(i));
  else
    fprintf('[%.2f,%.2f]', lb(i), ub(i));
  end
  if i < n
    if exist('fid', 'var')
      fprintf(fid, 'x');
    else
      fprintf('x');
    end
  end
end
end

