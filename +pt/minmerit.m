% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [it1, it2] = minmerit(it1, it2, opts)
% If it2 is empty, computes the merit function value for the current 
% iteration (it1.Fm) and a measure of its estimated decrease (it1.dm).
% If not, computes the merit function for the next iteration (it2.Fm).
%
% Assumes that all function values it1.fx, it1.ax, it1.aeqx, it.cx,
% it1.ceqx are already computed.
% Assumes that the square norm of the equalities it1.dceqx is already computed.
% Assumes that it1.d (a measure of the objectives decrease) is already
% computed.

% sizes
[~, nobj, ~, naeq, ~, nceq] = pt.mopsizes(it1);

if isempty(it2) % computations for it1
  if naeq + nceq == 0 % unconstrained case  
    it1.Fm = it1.fx; % the merit function is the objective function
    it1.dm = it1.d * ones(1, nobj);
    return
  end
    
  % equality constrained case
  if -it1.d <= opts.FirstOrdOptTol % can only make a progress relative to the equalities
    it1.Fm = 0.5 * it1.dceqx;
    it1.dm = -it1.dceqx;
  else 
    if it1.dceqx <= opts.FirstOrdOptTol % can only make a progress relative to the objectives
      it1.Fm = it1.fx; 
      it1.dm = it1.d * ones(1, nobj);
    else % can make a progress relative to the objectives and the equalities
      it1.Fm = [it1.fx, 0.5 * it1.dceqx];
      it1.dm = [it1.d * ones(1, nobj), -it1.dceqx];
    end
  end
else % computations for it2
  if naeq + nceq == 0 % unconstrained case  
    it2.Fm = it2.fx; % the merit function is the objective function
    return
  end
    
  % equality constrained case
  if -it1.d <= opts.FirstOrdOptTol % can only make a progress relative to the equalities
    it2.Fm = 0.5 * it2.dceqx;
  else 
    if it1.dceqx <= opts.FirstOrdOptTol % can only make a progress relative to the objectives
      it2.Fm = it2.fx; 
    else % can make a progress relative to the objectives and the equalities
      it2.Fm = [it2.fx, 0.5 * it2.dceqx];
    end
  end 
end
end