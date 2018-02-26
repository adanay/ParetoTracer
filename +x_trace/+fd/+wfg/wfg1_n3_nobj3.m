clear

% number of variables
n = 3;
nobj = 3;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.wfg.wfg1({n, nobj});

% real Pareto set and front
m = 10000;
[ps, m] = pset(m);
pf = objfun(ps);

% dimensions to be plotted
opts.PlotPSDims = 1 : 3;
opts.PlotPFDims = 1 : 3;

% plotting real Pareto set and front
utils.plotpareto(ps, pf, opts, 'linestyle', 'none', 'marker', '.', 'color', [0.8500, 0.3250, 0.0980]);
drawnow

% Calling the Pareto Tracer (PT) continuation algorithm.
x0 = ps(randi(m), :);

opts.HessApprox = 'fd';
opts.PCStepObj = 0.01;
opts.PCMinRelStepObj = -Inf;
opts.PCMaxRelStepObj = Inf;
opts.PCPlotMode = 'iter';

tic
[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
toc

