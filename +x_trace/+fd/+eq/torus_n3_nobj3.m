clear

% number of variables
n = 3;
nobj = 3;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.eq.torus({n, nobj});

% dimensions to be plotted
opts.PlotPSDims = 1 : 3;
opts.PlotPFDims = 1 : 3;

% Calling the Pareto Tracer (PT) continuation algorithm.
x0 = rand(1, n);
x0 = [0.1; -0.2; 0.4];

opts.HessApprox = 'fd';
objfun.H = [];
opts.PCStepObj = 0.15;
opts.LbObj = 10 * -ones(nobj, 1);
opts.UbObj = 10 * ones(nobj, 1);

tic
[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
toc

