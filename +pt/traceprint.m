% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = traceprint(info)

if strcmpi(info.opts.PCPrintMode, 'off')
  return
end

switch (upper(info.PHASEFLAG))
  case 'INIT'
    traceprintinit(info);
  case 'FINIT'
    traceprintfinit(info);
  otherwise
    traceprintit(info);
end

end

function [] = traceprintinit(info)

name = info.opts.MOPName;
if isempty(name)
  name = 'MOP';
end

[n, nobj, na, naeq, nc, nceq] = val.valmopsizes(info.sizes, false);

fprintf('\n');
fprintf([info.opts.PCIndent 'Pareto Tracer'], name, n, nobj);
fprintf('\n');
fprintf([info.opts.PCIndent 'Func: %s(n=%d,nobj=%d'], name, n, nobj);
if na > 0
  fprintf(',na=%d', na);
end
if naeq > 0
  fprintf(',naeq=%d', naeq);
end
if nc > 0
  fprintf(',nc=%d', nc);
end
if nceq > 0
  fprintf(',nceq=%d', nceq);
end
fprintf(')\n');

fprintf([info.opts.PCIndent 'Initial Point: [%s]\n'], utils.vector2str(info.it1.x, '%.2f'));
fprintf([info.opts.PCIndent 'Initial Fun Value: [%s]\n'], utils.vector2str(info.it1.fx, '%.2f'));

fprintf([info.opts.PCIndent 'Hess: Modif: %s, Approx: %s'], info.opts.HessModif, info.opts.HessApprox);
if strcmpi(info.opts.HessApprox, 'fd')
  if info.opts.FDForceHess
    full = 'Yes';
  else
    full = 'No';
  end
  fprintf([info.opts.PCIndent ', Approx Full: %s'], full);
end
fprintf([info.opts.PCIndent ' (Note: Approx used only if hess not provided.)']);
fprintf('\n');
if info.opts.PCUseStepVar
  fprintf([info.opts.PCIndent 'Step in Var: %.2f\n'], info.opts.PCStepVar);
else
  fprintf([info.opts.PCIndent 'Step in Obj: %.2f\n'], info.opts.PCStepObj);
end
end

function [] = traceprintfinit(info)

fprintf('\n');

% Exit
printtraceexitflag(info);

% Iterations
fprintf([info.opts.PCIndent 'Iterations: %d\n'], info.stats.PCIts);
fprintf([info.opts.PCIndent 'Solution Points: %d\n'], info.stats.Count);

%fprintf([info.opts.PCIndent 'Correct Stats: Avg Its: %.2f, Avg Dir Subprob Its: %.2f, Avg Lin Search Its: %.2f\n'], info.stats.OptIts, info.stats.OptDirIts, info.stats.OptLsIts);
fprintf([info.opts.PCIndent 'Correct Stats: Avg Its: %.2f, Avg Lin Search Its: %.2f\n'], info.stats.OptIts, info.stats.OptLsIts);

% Evals
pt.printevals(info, info.it1, info.opts.PCIndent);
end

function [] = traceprintit(info)

if ~strcmpi(info.opts.PCPrintMode, 'iter')
  return
end

switch(upper(info.PHASEFLAG))
  case 'IT-INIT'
    traceprintit_init(info);
  case 'PRED'
  case 'CORRECT'
  case 'IT-FINIT'
    traceprintit_finit(info);
end
end

function traceprintit_init(info)
fprintf([info.opts.PCIndent '%d: \n'], info.stats.PCIts);

% objective values
pt.printfunvals(info, info.it1, [info.opts.PCIndent info.opts.IndentGrowFactor]);
end

function traceprintit_finit(info)
indent = info.opts.IndentGrowFactor;

npred = length(info.itp);
ncorrect = length(info.itc);

% Predictors
if npred > 0
  predsteps = zeros(1, npred);
  for i = 1 : npred
    predsteps(i) = norm(info.itp(i).x - info.it1.x);
  end
  fprintf([info.opts.PCIndent indent 'Pred Steps: |f(p)-f(x)|=[%s]\n'], utils.vector2str(predsteps, '%.2f'));
end

% Correctors
if ncorrect > 0
  correctsteps = zeros(1, ncorrect);
  for i = 1 : ncorrect
    correctsteps(i) = norm(info.itc(i).x - info.it1.x);
  end
  fprintf([info.opts.PCIndent indent 'Correct Steps: |f(c)-f(x)|=[%s]\n'], utils.vector2str(correctsteps, '%.2f'));
end

%fprintf([info.opts.PCIndent indent 'Correct Stats: Avg Its: %.2f, Avg Dir Subprob Its: %.2f, Avg Lin Search Its: %.2f\n'], info.stats.LastOptIts, info.stats.LastOptDirIts, info.stats.LastOptLsIts);
fprintf([info.opts.PCIndent indent 'Correct Stats: Avg Its: %.2f, Avg Lin Search Its: %.2f\n'], info.stats.LastOptIts, info.stats.LastOptLsIts);

% KKT multipliers
pt.printkkt(info, info.it1, [info.opts.PCIndent indent]);

% Condition numbers
if pt.traceissd(info.objfun, info.multfun, info.opts)
  if ~info.it1.vIsSecant
    fprintf([info.opts.PCIndent indent 'rank/cond: W: %e\n'], info.it1.WRcond); % W = [J * Mu; [1...1 0...0]]
  end
else
  fprintf([info.opts.PCIndent indent 'rank/cond: Hw: %e, W: %e\n'], info.it1.HwxRcond, info.it1.WRcond); % W = [J * Mu; [1...1 0...0]]
end
end

function printtraceexitflag(info)
switch (info.EXITFLAG)
  case 0
    fprintf([info.opts.PCIndent 'Exit Flag (0): Number of iterations exceeded %d.\n'], info.opts.PCMaxIts);
  case 1
    fprintf([info.opts.PCIndent 'Exit Flag (1): No more solution points found.\n']);
  case -1
    fprintf([info.opts.PCIndent 'Exit Flag (-1): Stopped by an output function.\n']);
end
end


