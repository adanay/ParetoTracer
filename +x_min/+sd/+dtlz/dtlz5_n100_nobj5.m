clear

% number of variables
n = 100;
nobj = 5;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.dtlz.dtlz5({n, nobj});

% plotting real Pareto set and front
m = 10000;
ps = pset(m);
pf = objfun(ps);
utils.plotpareto(ps, pf, opts, 'linestyle', 'none', 'marker', '.', 'color', [0.8500, 0.3250, 0.0980]);
drawnow

sep = repmat('-', 1, 100); 

opts.HessApprox = 'off';

% Calling pt.minimize algorithm.
for i = 1 : 10
  x0 = rand(1, n);
  tic
  pt.minimize(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
  toc
  fprintf([sep '\n']);
end
