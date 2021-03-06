function [it] = minit()
% Initializes the structure that represents one iteration of a minimization
% algorithm.

% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% Lagrange multipliers
w = struct(...
  'objectives', [],...
  'lower', [],...
  'upper', [],...
  'ineqlin', [],...
  'eqlin', [],...
  'ineqnonlin', [],...
  'eqnonlin', []);
    
% previous iteration values
it = struct(...
  'x', [],...
  'fx', [],...
  'wfx', [],...
  'Jx', [],...
  'wJx', [],...
  'Hx', [],...
  'Hident', [],...
  'Lx', [],...
  'Lmodif', [],...
  'Hwx', [],...
  'ax', [],...
  'wax', [],...
  'aeqx', [],...
  'waeqx', [],...
  'cx', [],...
  'wcx', [],...
  'Jcx', [],...
  'wJcx', [],...
  'Hcx', [],...
  'Hcwx', [],...
  'ceqx', [],...
  'wceqx', [],...
  'Jceqx', [],...
  'wJceqx', [],...
  'Hceqx', [],...
  'Hceqwx', [],...
  'dceqx', [],... % square norm of the equalities
  'dcx', [],... % square norm of the inequalities
  'w', w,... % Lagrange multipliers
  'xActive', [],... % active dimensions regarding the box constraints
  'r', [],... % number of active dimensions regarding the box constraints
  'xLbActive', [],... % active dimensions regarding the lower bounds constraints
  'rLb', [],... % number of active dimensions regarding the lower bounds constraints
  'xUbActive', [],... % active dimensions regarding the upper bounds constraints
  'rUb', [],... % number of active dimensions regarding the upper bounds constraints
  'axActive', [],... % active linear inequalities
  'cxActive', [],... % active nonlinear inequalities 
  'v', [],... % search direction
  'd', [],... % measure of the objectives decrease
  'Fm', [],... % merit function
  'dm', [],... % measure of the merit function decrease
  't', [],... % step length
  'FirstOrdOpt', []);
end

