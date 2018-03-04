function [S, m] = usimplex(n, m)
% Uniform (n - 1)-simplex of approximately m points uniformly distributed.
% Each point component is between 0 and 1, i.e., 0 <= x1 <= 1, ..., 0 <= xn
% <= 1.
% Each point is a convex weight, i.e., x1 + x2 + ... + xn = 1.
% Returns the exact m required.

n1 = 1;
while nchoosek(n1 + n, n - 1) <= m
  n1 = n1 + 1;
end

S = nchoosek(1 : n1 + n - 1, n - 1) - repmat(0 : n - 2, nchoosek(n1 + n - 1, n - 1), 1) - 1;         
S = ([S, zeros(size(S, 1), 1) + n1] - [zeros(size(S, 1), 1), S]) / n1;

if n1 < n
  n2 = 0;
  while nchoosek(n1 + n - 1, n - 1) + nchoosek(n2 + n, n - 1) <= m
    n2 = n2 + 1;
  end
  if n2 > 0
    Sm = nchoosek(1 : n2 + n - 1, n - 1) - repmat(0 : n - 2, nchoosek(n2 + n - 1, n - 1), 1) - 1;
    Sm = ([Sm, zeros(size(Sm, 1), 1) + n2] - [zeros(size(Sm, 1), 1), Sm]) / n2;
    S  = [S; Sm / 2 + 1 / (2 * n)];
  end
end

S = max(S, 1e-6);
m = size(S, 1);
end

