function [it1, it2, stats, lsargout] = minls(it0, it1, objfun, lb, ub, lincon, nonlcon, multfun, opts, stats)
% Linear search using the Armijo condition with bisection.
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.
% Assumes that the square norm of the equalities it1.dceqx is already computed.
% Assumes that it1.d (a measure of the objectives decrease) is already
% computed.
% EXITFLAG: Describes the exit condition.
%  0: Number of iterations exceeded opts.OptLsMaxIts.
%  1: Armijo condition was satisfied.
% -2: No feasible step length was found.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% sizes
[~, ~, na, naeq, nc, nceq] = pt.mopsizes(it1);

[tl, tu] = lsinit();
zoom();

  function [] = zoom()
    its = 0;
    while tl == 0
      % no more progress is possible
      if tu <= sqrt(realmin)
        EXITFLAG = -2;
        it1.t = 0;
        it2 = it1;
        break
      end
      
      % max number of iterations reached
      if its == opts.OptLsMaxIts
        EXITFLAG = 0;
        % the new iteration cannot worsen the current solution
        if ~utils.isdominant(it2.Fm, it1.Fm)
          it1.t = 0;
          it2 = it1;
        end
        break;
      end
      
      nextit();
      its = its + 1;
      stats.OptLsIts = stats.OptLsIts + 1;
      
      % The Armijo condition is not satisfied on the merit function.
      if any(it2.Fm > it1.Fm + opts.ArmijoC * it1.t * it1.dm) || ~utils.isdominant(it2.Fm, it1.Fm)
        tu = it1.t;
        it1.t = 0.5 * (tl + tu); % bisection
      else
        EXITFLAG = 1;
        tl = it1.t;
      end
    end 
    lsargout = {its, EXITFLAG};
  end

  function [] = nextit()
    % Next iteration point.
    
    % next point of the line search
    it2.x = it1.x + it1.t * it1.v;
    
    % project onto the box
    it2.x = utils.project(it2.x, lb, ub);
    
    % function values
    [it2, stats] = pt.feval(it2, objfun, lincon, nonlcon, true, opts, stats);
    
    % computes the merit function value it2.Fm
    [~, it2] = pt.minmerit(it1, it2, opts);
  end

  function [tl, tu] = lsinit()
    % Initializes the linear search.
    
    % computes the merit function value it1.Fm and the estimated decrease it1.dm
    it1 = pt.minmerit(it1, [], opts);
    
    % the next iteration point
    it2 = pt.minit();
    
    % lower bound for t that satisfies the Armijo condition and the box constraints
    tl = 0; % since x must be feasible
    
    % upper bound for t that satisfies the box constraints
    s = zeros(size(it1.x));
    s = s + (it1.v < 0) .* (it1.x - lb);
    s = s + (it1.v > 0) .* (ub - it1.x);
    tu = min(s(s ~= 0) ./ abs(it1.v(s ~= 0)));
    
    % starting step length
    if pt.minissd(objfun, multfun, opts) && na + naeq + nc + nceq == 0 && stats.OptIts > 1
      t = it0.t * it0.d / it1.d;
      t = min(1, 1.01 * t);
      if t < eps^2
        t = 1;
      end
    else
      t = 1;
    end
    
    it1.t = min(t, tu);
  end
end


