function [isvalid] = valptoptnv(name, value, n)
% Validates the Pareto Tracer options.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if isempty(value)
  isvalid = true;
  return
end

switch(name)
  case 'HessApprox'
    isvalid = any(strcmpi({'off', 'bfgs', 'fd'}, value));
  case 'HessModif'
    isvalid = any(strcmpi({'off', 'chol'}, value));
  case 'FDForceHess'
    isvalid = isa(value, 'logical');
  case 'PCMaxIts'
    isvalid = isa(value, 'double') && value > 0;
  case 'OptMaxIts'
    isvalid = isa(value, 'double') && value > 0;
  case 'OptDirMaxIts'
    isvalid = isa(value, 'double') && value > 0;
  case 'OptDirSolver'
    isvalid = ischar(value);
  case 'OptLsMaxIts'
    isvalid = isa(value, 'double') && value > 0;
  case 'PCGMaxIts'
    isvalid = isa(value, 'double') && value > 0;
  case 'ArmijoC'
    isvalid = isa(value, 'double') && value > 0 && value < 1;
  case 'FirstOrdOptTol'
    isvalid = isa(value, 'double') && value > 0;
  case 'ConstViolTol'
    isvalid = isa(value, 'double') && value > 0;
  case 'OptMinStepVar'
    isvalid = isa(value, 'double') && value >= 0;
  case 'PCEdgeTol'
    isvalid = isa(value, 'double') && value > 0;
  case 'PCForceCells'
    isvalid = isa(value, 'logical');
  case 'PCStepObj'
    isvalid = isa(value, 'double');
  case 'PCMinRelStepObj'
    isvalid = isa(value, 'double');
  case 'PCMaxRelStepObj'
    isvalid = isa(value, 'double');
  case 'LbObj'
    isvalid = isa(value, 'double');
  case 'UbObj'
    isvalid = isa(value, 'double');
  case 'PCStepVar'
    isvalid = isa(value, 'double');
  case 'PCMinRelStepVar'
    isvalid = isa(value, 'double');
  case 'PCMaxRelStepVar'
    isvalid = isa(value, 'double');
  case 'PCUseStepVar'
    isvalid = isa(value, 'logical');
  case 'PCForceSecant'
    isvalid = isa(value, 'logical');
  case 'InitSetSize'
    isvalid = isa(value, 'double') && value > 0;
  case 'SetSizeGrowFact'
    isvalid = isa(value, 'double') && value > 1;
  case 'MOPName'
    isvalid = ischar(value);
  case 'OptOutFcn'
    isvalid = isa(value, 'function_handle') || ischar(value);
  case 'OptIndent'
    isvalid = ischar(value);
  case 'PCOutFcn'
    isvalid = isa(value, 'function_handle') || ischar(value);
  case 'PCIndent'
    isvalid = ischar(value);
  case 'IndentGrowFactor'
    isvalid = ischar(value);
  case 'SuppressOutput'
    isvalid = isa(value, 'logical');
  case 'PCSuppressOutput'
    isvalid = isa(value, 'logical');
  case 'PCSuppressPrint'
    isvalid = isa(value, 'logical');
  case 'PCSuppressPlot'
    isvalid = isa(value, 'logical');
  case 'OptSuppressOutput'
    isvalid = isa(value, 'logical');
  case 'OptSuppressPrint'
    isvalid = isa(value, 'logical');
  case 'OptSuppressPlot'
    isvalid = isa(value, 'logical');
  case 'OptPrintMode'
    isvalid = any(strcmpi({'off', 'iter', 'result'}, value));
  case 'PCPrintMode'
    isvalid = any(strcmpi({'off', 'iter', 'result'}, value));
  case 'OptPlotMode'
    isvalid = any(strcmpi({'off', 'iter', 'result'}, value));
  case 'PCPlotMode'
    isvalid = any(strcmpi({'off', 'iter', 'result', 'flow'}, value));
  case 'PCMarker'
    isvalid = true;
  case 'PCDrawCells'
    isvalid = isa(value, 'logical');
  case 'ValidateInput'
    isvalid = isa(value, 'logical');
  case 'PlotAxesPartition1'
    isvalid = isa(value, 'double') && value > 0;
  case 'PlotAxesPartition2'
    isvalid = isa(value, 'double') && value > 0;
  case 'PlotPSAxesInd'
    isvalid = isa(value, 'double') && value >= 0;
  case 'PlotPFAxesInd'
    isvalid = isa(value, 'double') && value >= 0;
  case 'PlotPSDims'
    isvalid = isa(value, 'double'); %&& length(value) <= 3 && all(value(:)' > 0, 2) && all(value(:)' <= n, 2);
  case 'PlotPFDims'
    isvalid = isa(value, 'double'); %&& length(value) <= 3 && all(value(:)' > 0, 2) && all(value(:)' <= n, 2);
  otherwise
    isvalid = val.valfdoptnv(name, value, n);
end
end
