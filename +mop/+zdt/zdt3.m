function [objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = zdt3(sizes, doval)
% ZDT3 multiobjective optimization problem.
% Suggested nobj = 2, n = 30.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% problem sizes
n = val.valmopsizes(sizes, true);
if n < 2
  n = 2;
end
nobj = 2;
sizes = struct(...
'variables', n,... 
'objectives', nobj);

if nargin < 2
  doval = false;
end

f1 = @(x) x(:, 1);
g = @(x) 1 + 9 * sum(x(:, 2 : n), 2) / (n - 1);

a = @(f1x, gx) 1 - sqrt(f1x ./ gx);
b = @(f1x) f1x .* sin(10 * pi * f1x);
c = @(f1x, gx) b(f1x) ./ gx;

h = @(f1x, gx) a(f1x, gx) - c(f1x, gx);

[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts] = mop.zdt.zdt(sizes, doval, f1, g, h);
pset = @zdt3pset;
pfront = @zdt3pfront;
opts.MOPName = 'ZDT3';

  function [ps, m] = zdt3pset(m)
    % disconnected segments of the ps
    lb = [0.0000000000, 0.1822287280, 0.4093136748, 0.6183967944, 0.8233317983];   
    ub = [0.0830015349, 0.2577623634, 0.4538821041, 0.6525117038, 0.8518328654];
    
    m5 = m / 5;
    
    ps = zeros(m, n); % each row is an individual
    ps(:, 1) = [linspace(lb(1), ub(1), m5),...
                linspace(lb(2), ub(2), m5),...
                linspace(lb(3), ub(3), m5),...
                linspace(lb(4), ub(4), m5),...
                linspace(lb(5), ub(5), m - 4 * m5)]; % for the ps set all variables to zero except the first one. 
  end

  function [pf, m] = zdt3pfront(m)
    % disconnected segments of the ps
    lb = [0.0000000000, 0.1822287280, 0.4093136748, 0.6183967944, 0.8233317983];   
    ub = [0.0830015349, 0.2577623634, 0.4538821041, 0.6525117038, 0.8518328654];
    
    m5 = m / 5;
    
    pf = [linspace(lb(1), ub(1), m5),...
          linspace(lb(2), ub(2), m5),...
          linspace(lb(3), ub(3), m5),...
          linspace(lb(4), ub(4), m5),...
          linspace(lb(5), ub(5), m - 4 * m5)]'; 
              
    pf = [pf, 1 - pf.^0.5 - pf .* sin(10 * pi * pf)];
  end
end
