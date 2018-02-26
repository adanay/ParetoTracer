% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [it1] = active(it1, lb, ub, force, opts)
% Determines the active inequality constraints.
% If force is true, the active constraints will be re-computed even if they
% are not empty.
% Assumes that all function values it1.fx, it1.ax, it.cx, are already computed.

if force || isempty(it1.xActive) || isempty(it1.xLbActive) || isempty(it1.xUbActive)
  [it1.xActive, it1.xLbActive, it1.xUbActive] = utils.isactive(it1.x, lb, ub, opts.ConstActiveTol);
  it1.r = sum(it1.xActive); % number of active dimensions respect to the box constraints
  it1.rLb = sum(it1.xLbActive); % number of active dimensions respect to the lower bounds constraints
  it1.rUb = sum(it1.xUbActive); % number of active dimensions respect to the upper bounds constraints
end

if force || isempty(it1.axActive)
  if ~isempty(it1.ax)
    it1.axActive = abs(it1.ax) < opts.ConstActiveTol;
  else
    it1.axActive = false(1, 0);
  end
end

if force || isempty(it1.cxActive)
  if ~isempty(it1.cx)
    it1.cxActive = abs(it1.cx) < opts.ConstActiveTol;
  else
    it1.cxActive = false(1, 0);
  end
end
end

