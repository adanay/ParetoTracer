clear

% number of variables
n = 100;
nobj = 2;

% obtaining test function handles
[objfun, lb, ub, lincon, nonlcon, multfun, sizes, opts, pset, pfront] = mop.eq.quad2({n, nobj});

% dimensions to be plotted
opts.PlotPSDims = [1, 100];
opts.PlotPFDims = 1 : 2;

sep = repmat('-', 1, 100); 

opts.HessApprox = 'off';
objfun.H = [];

% Calling pt.minimize algorithm.
for i = 1 : 10
  x0 = 0.1 * rand(1, n);
  tic
  pt.minimize(objfun, x0, [], lb, ub, lincon, nonlcon, multfun, opts);
  toc
  fprintf([sep '\n']);
end


