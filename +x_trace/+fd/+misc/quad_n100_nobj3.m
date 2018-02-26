clear

% number of variables
n = 100;
nobj = 3;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.misc.quad({n, nobj});

% real Pareto set and front
m = 10000;
[ps, m] = pset(m);
pf = objfun.f(ps);

% dimensions to be plotted
opts.PlotPSDims = 1 : 3;
opts.PlotPFDims = 1 : 3;

% plotting real Pareto set and front
%utils.plotpareto(ps, pf, opts, 'linestyle', 'none', 'marker', '.', 'color', [0.8500, 0.3250, 0.0980]);
%drawnow

% Calling the Pareto Tracer (PT) continuation algorithm.
x0 = rand(1, n);
%x0 = 0.5 * ones(n, 1);

opts.HessApprox = 'fd';
objfun.H = [];
opts.PCStepObj = 10;
opts.LbObj = 800 * -ones(nobj, 1);
opts.UbObj = 800 * ones(nobj, 1);

tic
[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
toc

