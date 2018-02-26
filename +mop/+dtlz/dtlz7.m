% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = dtlz7(sizes, doval)
% DTLZ7 multiobjective optimization problem.
% n >= 2, nobj >= 2, n >= nobj.
% Suggested nobj = 3, n = nobj + 19 = 22.

% problem sizes
[n, nobj] = val.valmopsizes(sizes, true);
if n < 2
  n = 2;
end
if nobj < 2
  nobj = 2;
end
if nobj > n
  nobj = n;
end

sizes = struct(...
  'variables', n,...
  'objectives', nobj);

if nargin < 2
  doval = false;
end

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts] = mop.dtlz.dtlz(sizes, @dtlz7f);
pset = @dtlz7pset;
pfront = @dtlz7pfront;
opts.MOPName = 'DTLZ7';

  function [fx] = dtlz7f(x)
    x = val.valmopargin(x, [], [], n, nobj, doval);
   
    y = x(:, 1 : nobj - 1);
    z = x(:, nobj : n);
    
    gz = 1 + 9 * mean(z, 2); 
        
    fx = [y, (1 + gz) .* (nobj - sum(y .* (1 + sin(3 * pi * y)), 2) ./ (1 + gz))];
  end

  function [ps, m] = dtlz7pset(m)
    q = n - nobj + 1;
    [ps, m] = utils.ugrid(nobj - 1, m);
    
    interval = [0, 0.251412, 0.631627, 0.859401];
    median = (interval(2) - interval(1)) / (interval(4) - interval(3) + interval(2) - interval(1));
    
    ps(ps <= median) = ps(ps <= median) * (interval(2) - interval(1)) / median + interval(1);
    ps(ps > median) = (ps(ps > median) - median) * (interval(4) - interval(3)) / (1 - median) + interval(3);
            
    ps = [ps, zeros(m, q)];
  end

  function [pf, m] = dtlz7pfront(m)
    [pf, m] = utils.ugrid(nobj - 1, m);
    
    interval = [0, 0.251412, 0.631627, 0.859401];
    median = (interval(2) - interval(1)) / (interval(4) - interval(3) + interval(2) - interval(1));
    
    pf(pf <= median) = pf(pf <= median) * (interval(2) - interval(1)) / median + interval(1);
    pf(pf > median) = (pf(pf > median) - median) * (interval(4) - interval(3)) / (1 - median) + interval(3);
            
    pf = [pf, 2 * (nobj - sum(pf / 2 .* (1 + sin(3 * pi * pf)), 2))];
  end
end



