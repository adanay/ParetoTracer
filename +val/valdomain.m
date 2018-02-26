% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Q] = valdomain(Q, n, doval)
% Validates a domain set.

if ~doval
  % no validation required
  Q = Q(:, :);
  return
end

if ~isa(Q, 'double')
  error('val:valdomain:NonDoubleQ', 'The domain Q must have double elements.');
end

if isempty(Q)
  Q = zeros(2, 0);
end

Q = Q(:, :);
mQ = size(Q, 1);

if mQ ~= 2
  error('val:valdomain:QBadShape', 'The domain Q must be of size (2 x n) where n is the number of variables.');
end

lb = Q(1, :);
ub = Q(2, :);
nQ = length(lb);
limit = 1e16;

if nQ > n
  warning('val:valdomain:IgnoringExtras', 'The domain Q has %d additional component(s) being ignored.', nQ - n);
  lb = lb(1 : n);
  ub = ub(1 : n);
elseif nQ < n
  if nQ > 0
    warning('val:valdomain:PadWithLimit', 'The domain Q has %d missing component(s) padded with (+-)%e.', n - nQ, limit);
  end
  lb = [lb, -limit * ones(1, n - nQ)];
  ub = [ub, limit * ones(1, n - nQ)];
end

if any(eq(ub, -inf))
  error('val:valdomain:MinusInfUb', 'The upper bounds of Q contains -Inf');
elseif any(eq(lb, inf))
  error('val:valdomain:PlusInfLb', 'The lower bounds of Q contains Inf.');
end

if any(isnan(lb))
  error('val:valdomain:NaNQ', 'The domain Q contains NaN.');
end
if any(isnan(ub))
  error('val:valdomain:NaNQ', 'The domain Q contains NaN.');
end

if ~isreal(lb)
  error('val:valdomain:NonRealQ', 'The domain Q contains nonreal entries.');
end
if ~isreal(ub)
  error('val:valdomain:NonRealQ', 'The domain Q contains nonreal entries.');
end
  
if any(lb > ub)
  error('val:valdomain:InfeasQ', 'The domain Q constraints provide an empty feasible set.');
end

lb(lb == -inf) = -limit;
ub(ub == inf) = limit;
Q = [lb; ub];
end


