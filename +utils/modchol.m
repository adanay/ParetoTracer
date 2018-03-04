function [L, modif] = modchol(A)
% Modified Cholesky decomposition to preserve positive definitness of A.
% Algorithm of Gill, Murray and Wright.
% Returns the lower triangular matrix L and whether or not it was modified.

tol = 5 * sqrt(eps);

n = length(A); % size of matrix
d = zeros(n, 1);
L = eye(n);
modif = 0;

diagA = diag(A);

gamma = max(abs(diagA));
xi = max(max(abs(A - diag(diagA))));

delta = eps * max(gamma + xi, 1);
beta = sqrt(max([gamma, xi / n, tol]));

for i = 1 : n
  p1 = 1 : i - 1;
  di = A(i, i) - L(i, p1) * (d(p1, 1) .* L(i, p1)');
  
  if i < n
    p2 = i + 1 : n;
    
    c = A(p2, i) - L(p2, p1) * (d(p1, 1) .* L(i, p1)');
    theta = max(abs(c));
    
    d(i) = max([abs(di), (theta / beta)^2, delta]);
    L(p2, i) = c / d(i);
  else
    d(i) = max(abs(di), delta);
  end
  
  if d(i) > di  % A was not sufficiently positive definite
    modif = 1;
  end
end
for i = 1 : n
  L(:, i) = L(:, i) * sqrt(d(i));
end
end

