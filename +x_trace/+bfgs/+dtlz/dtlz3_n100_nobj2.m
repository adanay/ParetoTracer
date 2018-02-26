clear

% number of variables
n = 100;
nobj = 2;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.dtlz.dtlz3({n, nobj});

% real Pareto set and front
m = 100;
[ps, m] = pset(m);
pf = pfront(m);

% dimensions to be plotted
opts.PlotPSDims = 1 : 3;
opts.PlotPFDims = 1 : 2;

% plotting real Pareto set and front
utils.plotpareto(ps, pf, opts, 'linewidth', 2, 'color', [0.8500, 0.3250, 0.0980]);
drawnow

% Calling the Pareto Tracer (PT) continuation algorithm.
x0 = ps(randi(m), :);

opts.PCStepObj = 0.04;

tic
[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
toc

