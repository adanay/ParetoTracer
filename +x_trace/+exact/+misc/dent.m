clear

n = 2;
nobj = 2;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.misc.dent();

% real Pareto set and front
m = 100;
[ps, m] = pset(m);
pf = objfun.f(ps);

% plotting real Pareto set and front
utils.plotpareto(ps, pf, opts, 'linewidth', 2, 'color', [0.8500, 0.3250, 0.0980]);
drawnow

% Calling the Pareto Tracer (PT) continuation algorithm.
x0 = rand(1, n);

opts.PCStepObj = 0.25;

tic
[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
toc

