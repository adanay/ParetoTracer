% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function[it1, stats] = mubase(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, force, opts, stats)
% Computes the mu vectors corresponding to an orthonormal basis of the 
% tangent space to the Pareto front at the current iteration point.
% It will compute (nobj - 1) vectors of size nobj.
% The components of mu will sum up to zero.
% If force is true, the vectors mu will be re-computed even if they are not 
% empty.

% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.
% Assumes that all objective and nonlinear constraint Jacobians it1.Jx, 
% it1.Jcx, it1.Jceqx are already computed.
% Assumes that the Lagrange multipliers w are already computed.

if ~force && ~isempty(it1.mu)
  return
end

% sizes
n = length(it1.x);
nobj = length(it1.fx);

% ensures the KKT multipliers are computed
[it1, stats] = pt.kkt(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, opts, stats);

% directions in objective space
% The tangent space to the Pareto front is orthogonal to the Lagrange
% multipliers corresponding to the objective functions.
[d, ~] = qr(it1.w.objectives'); 
d = d(:, 2 : nobj)';

% normalization
for i = 1 : nobj - 1
  normd = norm(d(i, :));
  if normd == 0
    continue
  end
  d(i, :) = d(i, :) / normd; 
end

% ensures the Mu space is computed
[it1, stats] = pt.Mu(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, opts, stats);

Jx = it1.Jx;
if it1.r < n % if r == n (corner) we will have to compute the unrestricted tangent and project later
  Jx = Jx(:, ~it1.xActive);
end

% computes mu based on d
[mu, R] = pt.mud(it1.Mu, Jx, d); 
it1.mu = mu;
it1.WRcond = R;
end