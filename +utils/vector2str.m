function [str] = vector2str(x, format)
% String representation of a vector.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if ~exist('format', 'var')
  format = '%e';
end

n = length(x); % number of variables
str = '';
for i = 1 : n
  str = strcat(str, sprintf(format, x(i)));
  if i < n
    str = strcat(str, ',');
  end
end
end

