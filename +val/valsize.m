function [s] = valsize(s)
% Validates an array of sizes by making it a vector of at least two
% elements where the last elements set to ones are trimmed.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

s = s(:)';
l = length(s);

switch(l)
  case 0
    s = [0 0];
    return
  case 1
    s = [s(1) 1];
    return
  case 2
    return
end

% number of ones at the end
count = 0;
for i = l : -1 : 3 % l must be 3 or more
  if s(i) == 1
    count = count + 1;
  else
    break
  end
end

% remove ones at the end
s = s(1 : l - count);
end

