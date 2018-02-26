clear

% number of variables
n = 100;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.zdt.zdt3({n});

% plotting real Pareto set and front
m = 1000;

ps1 = [0, 0, 0; 1, 0, 0];
pf1 = linspace(0, 1, m)';
pf1 = [pf1, 1 - pf1.^0.5 - pf1 .* sin(10 * pi * pf1)];
[ps2, m] = pset(m);
pf2 = pfront(m);

utils.plotpareto(ps1, pf1, opts, 'linewidth', 1, 'color', [0      0.4470 0.7410]);
utils.plotpareto(ps2, pf2, opts, 'linestyle', 'none', 'marker', '.', 'color', [0.8500 0.3250 0.0980]);
drawnow

sep = repmat('-', 1, 100);

% Calling pt.minimize algorithm.
for i = 1 : 10
  x0 = rand(1, n);
  tic
  pt.minimize(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
  toc
  fprintf([sep '\n']);
end