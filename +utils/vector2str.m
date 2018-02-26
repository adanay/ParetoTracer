% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [str] = vector2str(x, format)
% String representation of a vector.

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

