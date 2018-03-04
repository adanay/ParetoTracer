%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = minprint(info)

if strcmpi(info.opts.OptPrintMode, 'off')
  return
end

switch (upper(info.PHASEFLAG))
  case 'INIT' 
    minprintinit(info);
  case 'FINIT' 
    minprintfinit(info);
  otherwise
    minprintit(info);
end
end

function [] = minprintinit(info)

name = info.opts.MOPName;
if isempty(name)
  name = 'MOP';
end

[n, nobj, na, naeq, nc, nceq] = val.valmopsizes(info.sizes, false);

fprintf('\n');
fprintf([info.opts.OptIndent 'PT Minimize'], name, n, nobj);
fprintf('\n');
fprintf([info.opts.OptIndent 'Func: %s(n=%d,nobj=%d'], name, n, nobj);
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

fprintf([info.opts.OptIndent 'Initial Point: [%s]\n'], utils.vector2str(info.it1.x, '%.2f'));
fprintf([info.opts.OptIndent 'Initial Fun Value: [%s]\n'], utils.vector2str(info.it1.fx, '%.2f'));

fprintf([info.opts.OptIndent 'Hess: Modif: %s, Approx: %s'], info.opts.HessModif, info.opts.HessApprox);
if strcmpi(info.opts.HessApprox, 'fd')
  if info.opts.FDForceHess
    full = 'Yes';
  else
    full = 'No';
  end
  fprintf([info.opts.OptIndent ', Approx Full: %s'], full);
end
fprintf([info.opts.OptIndent ' (Note: Approx used only if hess not provided.)']);
fprintf('\n');
end

function [] = minprintfinit(info)

[~, ~, na, naeq, nc, nceq] = val.valmopsizes(info.sizes, false);

fprintf('\n');

% Exit
printminexitflag(info);

fprintf([info.opts.OptIndent 'Solution Point: [%s]\n'], utils.vector2str(info.result.x, '%.2f'));
fprintf([info.opts.OptIndent 'Solution Fun Value: [%s]\n'], utils.vector2str(info.result.fx, '%.2f'));

% Iterations
fprintf([info.opts.OptIndent 'Iterations: %d\n'], info.stats.OptIts);

% Optimality
fprintf([info.opts.OptIndent 'Optimality: d: %e, ||v||^2: %e'], info.result.d, info.result.v * info.result.v');
if na > 0 || nc > 0
  fprintf(', ||c||^2: %e', info.result.dcx);
end
if naeq > 0 || nceq > 0
  fprintf(', ||ceq||^2: %e', info.result.dcx);
end
fprintf('\n');

% DIR/LS Its
fprintf([info.opts.OptIndent 'Avg Dir Subprob Its: %.2f'], info.stats.OptDirIts);
fprintf(', Avg Lin Search Its: %.2f\n', info.stats.OptLsIts);

% Evals
pt.printevals(info, info.result, info.opts.OptIndent);
end

function [] = minprintit(info)

if ~strcmpi(info.opts.OptPrintMode, 'iter')
  return
end

indent = info.opts.IndentGrowFactor;

[~, ~, na, naeq, nc, nceq] = val.valmopsizes(info.sizes, false);

switch(upper(info.PHASEFLAG))
  case 'IT-INIT'
    % IT
    fprintf([info.opts.OptIndent '%d: \n'], info.stats.OptIts + 1);
    % objective values
    pt.printfunvals(info, info.it1, [info.opts.OptIndent info.opts.IndentGrowFactor]);
    
  case 'DIR'
    % DIR Its
    fprintf([info.opts.OptIndent indent 'Dir SubProb Its: %d\n'], info.OptDirIts);
    
    % Exit
    printdirexitflag(info);
    
    % Optimality
    fprintf([info.opts.OptIndent indent indent 'd: %e, ||v||^2: %e'], info.it1.d, info.it1.v * info.it1.v');
    if na > 0 || nc > 0
      fprintf(', ||c||^2: %e', info.it1.dcx);
    end
    if naeq > 0 || nceq > 0
      fprintf(', ||ceq||^2: %e', info.it1.dcx);
    end
    fprintf('\n');
    
    % KKT multipliers
    pt.printkkt(info, info.it1, [info.opts.OptIndent indent indent]);
    
  case 'STEP'
    % LS Its
    fprintf([info.opts.OptIndent indent 'Lin Search Its: %d\n'], info.OptLsIts);
 
    % Exit
    printlsexitflag(info);
    
    % Result
    fprintf([info.opts.OptIndent indent indent 't: %e\n'], info.it1.t);
end
end

function printminexitflag(info)
  [~, ~, na, naeq, nc, nceq] = val.valmopsizes(info.sizes, false);

  if na + naeq + nc + nceq > 0
    unconstrained = false;
  else
    unconstrained = true;
  end

  switch (info.EXITFLAG)
    case 0 
      fprintf([info.opts.OptIndent 'Exit (0): Number of iterations exceeded %d.\n'], info.opts.OptMaxIts);
    case 1 
      if unconstrained
        fprintf([info.opts.OptIndent 'Exit (1): First-order optimality measure was less than %e.\n'], info.opts.FirstOrdOptTol);
      else
        fprintf([info.opts.OptIndent 'Exit (1): First-order optimality measure was less than %e and constraint violation was less than %e.\n'], info.opts.FirstOrdOptTol, info.opts.ConstViolTol);
      end
    case 2 
      fprintf([info.opts.OptIndent 'Exit (2): Change in x too small (t = %e).\n'], info.result.t);
    case -1 
      fprintf([info.opts.OptIndent 'Exit (-1): Stopped by an output function.\n']);
    case -2 
      fprintf([info.opts.OptIndent 'Exit (-2): No feasible point was found.\n']);
  end
end

function printdirexitflag(info)
  indent = info.opts.IndentGrowFactor;
  
  switch(lower(info.opts.OptDirSolver))
    case 'fmincon'
      switch (info.DIREXITFLAG)
        case 0 
          fprintf([info.opts.OptIndent indent indent 'Exit (0): Number of iterations exceeded %d.\n'], info.opts.OptDirMaxIts);
        case 1 
          fprintf([info.opts.OptIndent indent indent 'Exit (1): First-order optimality measure was less than %e and constraint violation was less than %e.\n'], 1e-6, 1e-6);
        case 2 
          fprintf([info.opts.OptIndent indent indent 'Exit (2): Change in variables was less than %e and constraint violation was less than %e.\n'], 1e-6, 1e-6);
        case 3 
          fprintf([info.opts.OptIndent indent indent 'Exit (3): Change in objective function was less than %e and constraint violation was less than %e.\n'], 1e-6, 1e-6);
        case -1 
          fprintf([info.opts.OptIndent indent indent 'Exit (-1): Stopped by an output function.\n']);
        case -2 
          fprintf([info.opts.OptIndent indent indent 'Exit (-2): No feasible point was found.\n']);
        case -3
          fprintf([info.opts.OptIndent indent indent 'Exit (-3): Objective function went below %e and constraint violation was less than %e.\n'], -1.0000e+20, 1e-6);
        otherwise
          fprintf([info.opts.OptIndent indent indent 'Exit (%d): Unknown.\n'], info.DIREXITFLAG);
      end
  end
end

function printlsexitflag(info)
  indent = info.opts.IndentGrowFactor;
  
  switch (info.LSEXITFLAG)
    case 0 
      fprintf([info.opts.OptIndent indent indent 'Exit (0): Number of iterations exceeded %d.\n'], info.opts.OptLsMaxIts);
    case 1 
      fprintf([info.opts.OptIndent indent indent 'Exit (1): Armijo condition was satisfied.\n']);
    case -2 
      fprintf([info.opts.OptIndent indent indent 'Exit (-2): No feasible step length was found.\n']);
  end
end



