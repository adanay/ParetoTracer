clear

% number of variables
n = 100;
nobj = 2;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.eq.quad2({n, nobj});

% dimensions to be plotted
opts.PlotPSDims = [1, 100];
opts.PlotPFDims = 1 : 2;

% Calling the Pareto Tracer (PT) continuation algorithm.
x0 = 0.1 * rand(n, 1);

objfun.H = [];
opts.PCStepObj = 0.05;
opts.PCPlotMode = 'iter';

tic
[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
toc

