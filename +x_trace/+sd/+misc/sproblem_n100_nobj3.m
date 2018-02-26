clear

% number of variables
n = 100;
nobj = 3;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.misc.sproblem({n, nobj});

% dimensions to be plotted
opts.PlotPSDims = 1 : 3;
opts.PlotPFDims = 1 : 3;

% Calling the Pareto Tracer (PT) continuation algorithm.
x0 = rand(1, n);
%x0 = 0.5 * ones(n, 1);

opts.HessApprox = 'off';
objfun.H = [];
opts.PCStepObj = 10;
opts.LbObj = 800 * -ones(nobj, 1);
opts.UbObj = 800 * ones(nobj, 1);

tic
[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
toc

