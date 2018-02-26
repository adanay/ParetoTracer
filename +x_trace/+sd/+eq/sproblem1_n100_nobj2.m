clear

% number of variables
n = 100;
nobj = 2;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.eq.sproblem1({n, nobj});

% dimensions to be plotted
opts.PlotPSDims = 1 : 3;
opts.PlotPFDims = 1 : 2;

% Calling the Pareto Tracer (PT) continuation algorithm.
x0 = rand(1, n);

opts.HessApprox = 'off';
objfun.H = [];
opts.PCStepObj = 10;

tic
[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
toc

