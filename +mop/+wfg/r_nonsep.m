function [value] = r_nonsep(y, i, A)
% A in {1, 2,...,|y|}.
% |y| mod A = 0.
% i represents the considered indices from y.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

y = y(:, i);
n = size(y, 2);

value = 0;
for j = 1 : n
  value = value + y(:, j);
  for k = 0 : A - 2
    value = value + abs(y(:, j) - y(:, 1 + mod(j + k, n)));
  end
end

c = n * ceil(A / 2) * (1 + 2 * A - 2 * ceil(A / 2)) / A;
value = value / c;
end

