% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [it1, stats] = tangent(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, force, objToMin, opts, stats)
% Computes v, the tangent(s) to the Pareto set. v is of size (m x n) where
% n is the number of variables and m is the number of tangent vectors,  
% which corresponds to the number of vectors it.mu (if not empty). 
%
% v = -Mu * mu.
%
% If it.mu is empty, it will compute (nobj - 1) vectors corresponding to an 
% orthonormal basis of the tangent space to the Pareto front at the current 
% iteration point.  
%
% For the bi-objective case, it may end up computing secants instead of 
% tangents if no Hessian information is available. For this, it0.x must not
% be empty since the secant uses the previous iteration point. Additionally, 
% it1.mu must be empty (if it is provided, the method uses it to compute
% v = -Mu * mu). The secant vectors are normalized. If secants are used, 
% it1.vIsSecant is set to true such that further steps can be done accordingly. 
%
% For the bi-objective case as well, if objToMin = i, the tangent direction  
% will be adjusted to minimize the i-th objective.
%
% If force is true, the tangent space will be re-computed even if it is not 
% empty.
%
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.
% Assumes that all objective and nonlinear constraint Jacobians it1.Jx, 
% it1.Jcx, it1.Jceqx are already computed.
% Assumes that the Lagrange multipliers w are already computed.

if ~force && ~isempty(it1.v)
  return
end

% sizes
n = length(it1.x);
nobj = length(it1.fx);
if nobj == 1
  it1.v = zeros(0, n);
  return;
end

% active sets 
it1 = pt.active(it1, lb, ub, false, opts);

% Since a linear combination of derivatives may be nearly orthogonal to
% the PS, the secant seems to be a better choice than derivatives when 
% there is no Hessian information.
% Applies only to the bi-objective case (where it1.mu is empty).
if nobj == 2 && ~isempty(it0.x) &&...
  (isempty(it1.mu) && pt.traceissd(objfun, multfun, opts) ||...
   opts.PCForceSecant)
    
  % Excluding active box constraints since the tangent will have to be zero  
  % in those components.
  if it1.r < n % if r == n (corner) we will have to compute the unrestricted tangent and project later
    v = zeros(1, n); 
    v(~it1.xActive) = it1.x(~it1.xActive) - it0.x(~it1.xActive); % secant
  else
    v = it1.x - it0.x; % secant
  end
  
  % normalization
  normv = norm(v);
  if normv ~= 0
    v = v / normv;
    it1.vIsSecant = true;
  else
    tangent_Mu()
  end
else 
  tangent_Mu()
end

% number of tangent vectors
m = size(v, 1);

% direction adjustment
if nobj == 2 && m == 1 && ~isempty(objToMin) && objToMin ~= 0
  it1.Jvx = it1.Jx * v';
  stats.JvCount = stats.JvCount + 1;
  
  if it1.Jvx(objToMin) > 0
    v = -v;
  end 
end

it1.v = v;

  function tangent_Mu()
    if isempty(it1.mu)
      % computes (nobj - 1) mu vectors corresponding to an orthonormal
      % basis of the tangent space to the PF at x.
      [it1, stats] = pt.mubase(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, opts, stats);
    end
    
    % Ensures that the Mu space is computed.
    [it1, stats] = pt.Mu(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, false, opts, stats);
    
    % Excluding active box constraints since the tangent will have to be zero
    % in those components.
    if it1.r < n % if r == n (corner) we will have to compute the unrestricted tangent and project later
      v = zeros(size(it1.mu, 1), n);
      v(:, ~it1.xActive) = -(it1.Mu * it1.mu')';
    else
      v = -(it1.Mu * it1.mu')';
    end
  end
end


