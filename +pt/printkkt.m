%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = printkkt(info, it, indent)
[~, ~, na, naeq, nc, nceq] = pt.mopsizes(it);

w = it.w;

% KKT multipliers
fprintf([indent 'KKT: objs: [%s]=%.2f\n'], utils.vector2str(w.objectives, '%.2f'), sum(w.objectives));
if na > 0
  fprintf([indent '     ineqlin: [%s]\n'], utils.vector2str(w.ineqlin, '%.2f'));
end
if naeq > 0
  fprintf([indent '     eqlin: [%s]\n'], utils.vector2str(w.eqlin, '%.2f'));
end
if nc > 0
  fprintf([indent '     ineqnonlin: [%s]\n'], utils.vector2str(w.ineqnonlin, '%.2f'));
end
if nceq > 0
  fprintf([indent '     eqnonlin: [%s]\n'], utils.vector2str(w.eqnonlin, '%.2f'));
end
if it.rLb > 0
  fprintf([indent '     lower: [%s]\n'], utils.vector2str(w.lower, '%.2f'));
end
if it.rUb > 0
  fprintf([indent '     upper: [%s]\n'], utils.vector2str(w.upper, '%.2f'));
end
end

