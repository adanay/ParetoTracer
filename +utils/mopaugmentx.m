% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [x] = mopaugmentx(x, w, maxdim)
n = length(x);
nobj = length(w.objectives);
na = length(w.ineqlin);
naeq = length(w.eqlin);
nc = length(w.ineqnonlin);
nceq = length(w.eqnonlin);

if maxdim > n
  x = [x, w.objectives];
  if maxdim > n + nobj
    x = [x, w.ineqlin];
    if maxdim > n + nobj + na
      x = [x, w.eqlin];
      if maxdim > n + nobj + na + naeq
        x = [x, w.ineqnonlin];
        if maxdim > n + nobj + na + naeq + nc
          x = [x, w.eqnonlin];
          if maxdim > n + nobj + na + naeq + nc + nceq
            x = [x, w.lower];
            if maxdim > n + nobj + na + naeq + nc + nceq + n
              x = [x, w.upper];
            end
          end
        end
      end
    end
  end
end
end