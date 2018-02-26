% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = printfunvals(info, it, indent)
[~, ~, na, naeq, nc, nceq] = pt.mopsizes(it);

fprintf([indent 'f(%s)=[%s]\n'], utils.vector2str(it.x, '%.2f'), utils.vector2str(it.fx, '%.2f'));
if na > 0
  fprintf([indent 'A*x-b=[%s]\n'], utils.vector2str(it.ax, '%.2f'));
end
if naeq > 0
  fprintf([indent 'Aeq*x-beq=[%s]\n'], utils.vector2str(it.aeqx, '%.2f'));
end
if nc > 0
  fprintf([indent 'c(%s)=[%s]\n'], utils.vector2str(it.x, '%.2f'), utils.vector2str(it.cx, '%.2f'));
end
if nceq > 0
  fprintf([indent 'ceq(%s)=[%s]\n'], utils.vector2str(it.x, '%.2f'), utils.vector2str(it.ceqx, '%.2f'));
end
end

