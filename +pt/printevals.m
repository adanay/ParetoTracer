%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = printevals(info, it, indent)
[~, ~, na, naeq, nc, nceq] = pt.mopsizes(it);

fprintf('\n');

% Fun Evals
fprintf([indent 'Fun Evals: %d\n'], info.stats.fCount);

if na > 0
  fprintf([indent 'Lin Ineq Evals: %d'], info.stats.aCount);
end
if naeq > 0
  if na > 0
    fprintf(',  Lin Eq Evals: %d', info.stats.aeqCount);
  else
    fprintf([indent 'Lin Eq Evals: %d'], info.stats.aeqCount);
  end
end

if na + naeq > 0
  fprintf('\n');
end

if nc > 0
  fprintf([indent 'Nonlin Ineq Evals: %d'], info.stats.cCount);
end
if nceq > 0
  if nc > 0
    fprintf(',  Nonlin Eq Evals: %d', info.stats.ceqCount);
  else
    fprintf([indent 'Nonlin Eq Evals: %d'], info.stats.ceqCount);
  end
end

if nc + nceq > 0
  fprintf('\n');
end

% Jac Evals
if isstruct(info.objfun) && ~isempty(info.objfun.J)
  fprintf([indent 'Jac Evals: %d\n'], info.stats.JCount);
end

temp = false;
if ~isempty(info.nonlcon) && isstruct(info.nonlcon) && ~isempty(info.nonlcon.Jc)
  temp = true;
  fprintf([indent 'Jac Nonlin Ineq Evals: %d'], info.stats.JcCount);
end
if ~isempty(info.nonlcon) && isstruct(info.nonlcon) && ~isempty(info.nonlcon.Jceq)
  if temp == true
    fprintf(',  Jac Nonlin Eq Evals: %d', info.stats.JceqCount);
  else
    fprintf([indent 'Jac Nonlin Eq Evals: %d'], info.stats.JceqCount);
    temp = true;
  end
end
if temp == true
  fprintf('\n');
end

% Hess Evals
if isstruct(info.objfun) && ~isempty(info.objfun.H)
  fprintf([indent 'Hess Evals: %d\n'], info.stats.HCount);
end
end

