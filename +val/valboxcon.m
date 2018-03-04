function [lb, ub] = valboxcon(lb, ub, n, doval)
% Validates the box constraints of an optimization solver.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

limit = 1e16;

if ~doval
  % no validation required
  lb = lb(:)';
  ub = ub(:)';
  if isempty(lb)
    lb = -limit * ones(1, n);
  end
  if isempty(ub)
    ub = limit * ones(1, n);
  end
  return
end

if ~isa(lb, 'double') || ~isa(ub, 'double')
  error('val:valboxcon:NonDoubleLbUb', 'The vectors lb and ub must have double elements.');
end

lb = lb(:)';
ub = ub(:)';
nlb = length(lb);
nub = length(ub);

if nlb > n
  warning('val:valboxcon:IgnoringExtraLbs', 'The lower bounds lb has %d additional component(s) being ignored.', nlb - n);
  lb = lb(1 : n);
elseif nlb < n
  if nlb > 0
    warning('val:valboxcon:PadLbWithMinusLimit', 'The lower bounds lb has %d missing component(s) padded with -%e.', n - nlb, limit);
  end
  lb = [lb, -limit * ones(1, n - nlb)];
end

if nub > n
  warning('val:valboxcon:IgnoringExtraUbs', 'The upper bounds ub has %d additional component(s) being ignored.', nub - n);
  ub = ub(1 : n);
elseif nub < n
  if nub > 0
    warning('val:valboxcon:PadUbWithLimit', 'The upper bounds ub has %d missing component(s) padded with %e.', n - nub, limit);
  end
  ub = [ub, limit * ones(1, n - nub)];
end

if any(eq(ub, -inf))
  error('val:valboxcon:MinusInfUb', 'The upper bounds ub contains -Inf.');
elseif any(eq(lb, inf))
  error('val:valboxcon:PlusInfLb', 'The lower bounds lb contains Inf.');
end

if any(isnan(lb))
  error('val:valboxcon:NaNLb', 'The lower bounds lb contains NaN.');
end
if any(isnan(ub))
  error('val:valboxcon:NaNUb', 'The upper bounds ub contains NaN.');
end

if ~isreal(lb)
  error('val:valboxcon:NonRealLb', 'The lower bounds lb contains nonreal entries.');
end
if ~isreal(ub)
  error('val:valboxcon:NonRealUb', 'The upper bounds ub contains nonreal entries.');
end
  
if any(lb > ub)
  error('val:valboxcon:InfeasLbUb', 'The box constraints provide an empty feasible set.');
end

lb(lb == -inf) = -limit;
ub(ub == inf) = limit;
end

