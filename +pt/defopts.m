function [opts] = defopts(n, nobj)
% Default options for Pareto Tracer.
% - HessApprox: Ignored if the Hessian or the Hessian multiply functions 
%   are provided.
%   + 'bfgs' (default): A Quasi-Newton BFGS update is used for each Hessian. 
%   + 'off': A steepest descent method is performed, i.e., the Hessians are
%     assumed to be the identity matrix.
%   + 'fd': Finite differences are used to approximate all the Hessians.
% - HessModif: Defines the type of modification applied to the Hessians to 
%   ensure they are positive definite.
%   + 'chol': Modified Cholesky decomposition. This is the default.
%   + 'off': No modification is performed.
% - FDForceHess: When using FD, defines whether to approximate a full Hessian 
%   at each iteration instead of several products v' * H. By default it is 
%   false or true if n <= 5. If opts.LargeScale was specified, this setting 
%   has no effect as the idea is not to form (n x n) matrices.
% - FunValCheck: Checks whether the objective function values are valid. If 
%   true, displays an error when the objective function returns a value  
%   that is complex, NaN, or Inf. If the Jacobian is not valid, a warning 
%   will be thrown and the Jacobian will be approximated by FD. Same for
%   the Hessians, but they can also be updated by QN, or approximated to the 
%   identity (depending on the HessApprox choice). The default is true.
% - LargeScale: If the problem is large scale, a matrix of size (n x n)
%   will never be formed unless the Hessian is provided, the multiply 
%   function Hw is provided, or opts.HessApprox = 'bfgs'. Only matrices of 
%   size (nobj x n) (same size as the Jacobian) will be formed. 
%   The default is false.  
% - PCMaxIts: 10000. Max number of iterations allowed to the continuation  
%   algorithm, where each iteration consists of 
%   - a predictor stage (one or several predictors around a point),
%   - and a corrector (optimization) stage (those predictors are corrected).
% - PCGMaxIts: 20. Max no of iterations allowed to the PCG algorithm utilized 
%   to compute the Mu space. 
% - OptMaxIts: 1000. Max number of iterations allowed to optimization 
%   algorithms.
% - OptDirSolver: fmincon. Solver utilized to solve the optimization 
%   direction subproblem.
% - OptDirMaxIts: 20. Max number of iterations allowed to the optimization 
%   direction subproblem.
% - OptLsMaxIts: 20. Max number of iterations allowed to the linear search 
%   strategy of optimization algorithms.
% - ArmijoC: 0.1. Constant utilized by linear searches in optimization
%   algorithms.
% - FirstOrdOptTol: 1e-6. First order optimality tolerance.
% - OptMinStepVar: 1e-10. Lower bound on the length of a step. The
%   optimization algorithm will stop if the newly computed point is too
%   close to the previous one.
% - ConstViolTol: 1e-6. Tolerance on the constraint violation.
% - ConstActiveTol: 1e-6. Tolerance utilized to determine whether or not a 
%   variable is active with respect to the contraints.
% - PCEdgeTol: 1e-3. If there is a component of alpha below this tolerance,
%   the current iteration point is considered to be on an edge.
% - PCStepObj: Scalar step length in objective space. The default is 0.1.
% - PCMinRelStepObj: If the relative step length (regarding PCStepObj) is
%   below this tolerance, the new point will be discarded.
% - PCMaxRelStepObj: If the relative step length (regarding PCStepObj) is
%   above this tolerance, the new point will be discarded.
% - LbObj: Vector representing the lower bounds in objective space.
% - UbObj: Vector representing the upper bounds in objective space.
% - PCForceCells: Force the use of space cells even if the problem is 
%   bi-objective. Space cells are mandatory for nobj > 2. 
% - PCForceSecant: Force the use of secants for bi-objective problems.
% - InitSetSize: Initial array length for variable length arrays, e.g., the 
%   computed solution set of a PC algorithm.
% - SetSizeGrowFact: Grow factor for variable length arrays.
% - OptOutFcn: Function to display info or to stop the optimization  
%   algorithm. It has the following format:
%   function [stop, it1, it2, stats] = optoutfcn(info)
%   where info is a structure containing the current variables:
%   + it0, it1, it2: Previous, current, and next iteration structures. See
%     pt.minit for the list of fields.
%   + objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts: The
%     input parameters to the optimization function.
%   + result, stats, EXITFLAG: The output parameters of the optimization
%     function.
%   + PHASEFLAG = 'INIT': The algorithm is ready to start.
%               = 'IT-INIT': The current iteration is ready to start.
%               = 'DIR': The descent direction was just computed.
%               = 'STEP': The step length was just computed.
%               = 'IT-FINIT': The current iteration is about to finish.
%               = 'FINIT': The algorithm has finished.
%   + DIREXITFLAG: Exit condition of the direction subproblem algorithm.
%   + OptDirIts: No of iterations performed by the direction subproblem
%     solver.
%   + LSEXITFLAG: Exit condition of the line search algorithm.
%   + OptLsIts: No of iterations performed by the line search.
%   The returned values (if any) replace the values utilized by the algorithm.
%
% - PCOutFcn: Function to display info or to stop the PC algorithm. It is 
%   also used to discard the current predictor or corrector point. It has 
%   the following format:
%   function [stop, discard, it1, itp, itc, stats] = pcoutfcn(info)
%   where info is a structure containing the current variables:
%   + it0, it1: Previous and current iteration structures. See pt.traceit 
%     for the list of fields.
%   + itp, itc: Current predictor or corrector structures. See pt.traceit 
%     for the list of fields.  
%   + objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts: The
%     input parameters to the continuation function.
%   + result, stats, EXITFLAG: The output parameters of the continuation
%     function.
%   + PHASEFLAG = 'INIT': The algorithm is ready to start.
%               = 'IT-INIT': The current iteration is ready to start.
%               = 'PRED': A new predictor was just computed.
%               = 'CORRECT': A new corrector was just computed.
%               = 'SAVE': The corrector was saved.  
%               = 'IT-FINIT': The current iteration is about to finish.
%               = 'FINIT': The algorithm has finished.
%   The returned values (if any) replace the values utilized by the algorithm.
%
% The following options are utilized by output functions only.
%
% - MOPName: Name of the multiobjective optimization problem (MOP) being
%   solved. 
% - PCIndent, OptIndent: Initial indentation for the output of the algorithms. 
%   By default it is ''.
% - IndentGrowFactor: It is '  ' by default.
% - PCSuppressOutput, PCSuppressPrint, PCSuppressPlot, 
%   OptSuppressOutput, OptSuppressPrint, OptSuppressPlot: Controls the output.
%   They are all false by default.
% - OptPrintMode: 'off', 'iter', 'result'. By default is 'result'.
% - PCPrintMode: 'off', 'iter', 'result'. By default is 'result'.
% - OptPlotMode: 'off', 'iter', 'result'. By default is 'iter'.
% - PCPlotMode: 'off', 'iter', 'result', 'flow'. By default is 'result'.
%   For biobjective problems it is 'flow'.
% - PCDrawCells: true/false. Determines whether to draw the computed cells.
%   By default it is false. 
% - PlotAxesPartition1, PlotAxesPartition2: Breaks the Figure window into 
%   an m-by-n matrix of small axes. By default it is (1 x 2).
% - PlotPSAxesInd, PlotPFAxesInd: Indices of the axes where the PS and the  
%   PF will be plotted. 0 means that they are not plotted. By default they
%   are 2 and 1 respectively.
% - PlotPSDims, PlotPFDims: Respectively the dimensions of the PS and PF to
%   be plotted. By default they are both 1 : 3.
%
% Additionally, all the options for the computation of finite differences
% can be specified. See fd.defopts for the list of available choices.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

opts = struct(...
  'HessApprox', 'bfgs',...% (bfgs, fd, off)
  'HessModif', 'chol',...% (chol, off)
  'FDForceHess', false,...% Force full Hessian approximations instead of 
                       ...% directional derivatives (H * v).  
  'PCMaxIts', 10000,...% Max its done by pt.trace (pred-correct its).
  'PCGMaxIts', 20,...% Max its done by the pred dir solver(s) (pcg). 
  'OptMaxIts', 100,...% Max its done by pt.minimize.
  'OptDirSolver', 'fmincon',...
  'OptDirMaxIts', 20,...% Max its done by the opt dir solver (fmincon).
  'OptLsMaxIts', 20,...% Max its done by the line search.
  'ArmijoC', 0.1,...
  ...
  'FirstOrdOptTol', 1e-6,...% First order optimality tol.
  'OptMinStepVar', 1e-10,... % pt.minimize will stop if the step is below this tol.
  'ConstViolTol', 1e-6,...% Constraint violation tol.
  'ConstActiveTol', 1e-6,...% Tol for considering that a constraint is active.
  'PCEdgeTol', 1e-3,...% Tol for considering that a point is at an edge of the PS/PF.
  ...
  'PCStepObj', 1e-1,...% Step length in objective space (distance between computed optimal points).
  'PCMinRelStepObj', 1e-1,... % pt.trace will discard points where the relative step is below this tol.
  'PCMaxRelStepObj', 1e+1,... % pt.trace will discard points where the relative step is above this tol.
  'LbObj', [],...% Lower bounds in objective space.
  'UbObj', [],...% Upper bounds in objective space.
  ...
  'PCStepVar', 1e-1,...% Not implemented. Step length in decision space (distance between computed optimal points).
  'PCMinRelStepVar', 1e-1,... % Not implemented. pt.trace will discard points where the relative step is below this tol.
  'PCMaxRelStepVar', 1e+1,... % Not implemented. pt.trace will discard points where the relative step is above this tol.
  'PCUseStepVar', false,...% Not implemented. Use the step var given in decision space.   
  ...
  'PCForceSecant', false,...% Force the use of secants for bi-objective problems.
  'PCForceCells', false,...% Use cells for bi-objective problems.
  ...
  'InitSetSize', 100,...
  'SetSizeGrowFact', 1.5,...
  ...
  'OptOutFcn', @pt.minout,...
  'PCOutFcn', @pt.traceout,...
  ...% Below all options for the out functions.
  'SuppressOutput', false,...
  'PCSuppressOutput', false,...
  'PCSuppressPrint', false,...
  'PCSuppressPlot', false,...
  'OptSuppressOutput', false,...
  'OptSuppressPrint', false,...
  'OptSuppressPlot', false,...
  ...
  'MOPName', 'MOP',...
  'OptIndent', '',...
  'PCIndent', '',...
  'IndentGrowFactor', '  ',...
  'OptPrintMode', 'result',...% (result, iter, off)
  'PCPrintMode', 'result',...% (result, iter, off)
  ...
  'OptPlotMode', 'iter',...% (result, iter, off)
  'PCPlotMode', 'result',...% (result, flow, iter, off)
  'PCMarker', 'o',...
  'PCDrawCells', false,...% draws the computed cells 
  'PlotAxesPartition1', 1,...% horizontal partitions of the figure
  'PlotAxesPartition2', 2,...% vertical partitions of the figure
  'PlotPSAxesInd', 2,...% partition index where the PS will be plotted 
                     ...% (0 means do not plot)
  'PlotPFAxesInd', 1,...% partition index where the PF will be plotted
                     ...% (0 means do not plot)
  'PlotPSDims', 1 : 3,...% variable indices to be plotted (up to 3)
  'PlotPFDims', 1 : 3);  % objective indices to be plotted (up to 3) 

if nargin > 0
  if n <= 5
    opts.FDForceHess = true;
  end
  opts.PlotPSDims = 1 : min(3, n);
end

if nargin > 1
  if nobj > 2
    opts.PCMarker = 'o';
  else
    opts.PCPlotMode = 'flow';
  end
  opts.PlotPFDims = 1 : min(3, nobj);
end

fdopts = fd.defopts();
fdopts.UseVectorized = false;

names = [fieldnames(opts); fieldnames(fdopts)];
opts = cell2struct([struct2cell(opts); struct2cell(fdopts)], names, 1);
end

