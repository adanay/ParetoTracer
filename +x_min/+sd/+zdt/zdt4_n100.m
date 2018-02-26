clear

% number of variables
n = 100;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.zdt.zdt4({n});

% plotting real Pareto set and front
m = 100;
ps = pset(m);
pf = pfront(m);
utils.plotpareto(ps, pf, opts, 'linewidth', 2, 'color', [0.8500, 0.3250, 0.0980]);
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
