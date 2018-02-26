clear

% number of variables
n = 100;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.zdt.zdt3({n});

% real Pareto set and front
m = 1000;
[ps, m] = pset(m);
pf1 = linspace(0, 1, m)';
pf1 = [pf1, 1 - pf1.^0.5 - pf1 .* sin(10 * pi * pf1)];
pf2 = pfront(m);

% dimensions to be plotted
opts.PlotPSDims = 1 : 3;
opts.PlotPFDims = 1 : 2;

% plotting real Pareto set and front
utils.plotpareto(ps, pf1, opts, 'linewidth', 2, 'color', [0      0.4470 0.7410]);
utils.plotpareto(ps, pf2, opts, 'linestyle', 'none', 'marker', '.', 'color', [0.8500 0.3250 0.0980]);
drawnow

% Calling the Pareto Tracer (PT) continuation algorithm.
x0 = ps(randi(m), :);

opts.HessApprox = 'off';
opts.PCStepObj = 0.02;
opts.PCPlotMode = 'result';

tic
[result, stats, EXITFLAG] = pt.trace(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
toc

