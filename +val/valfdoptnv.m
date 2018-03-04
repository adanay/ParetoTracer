function [isvalid] = valfdoptnv(name, value, n)
% Validates the FD options.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if isempty(value)
  isvalid = true;
  return
end

switch(name)
  case 'FDType'
    isvalid = any(strcmpi({'forward', 'central'}, value));
  case 'TypicalX'
    isvalid = isa(value, 'double') && (length(value) == 1 || length(value) == n);
  case 'FDMinChange'
    isvalid = isa(value, 'double');
  case 'FDMaxChange'
    isvalid = isa(value, 'double');
  case 'FDStepSize'
    isvalid = isa(value, 'double') && (length(value) == 1 || length(value) == n);
  case 'FDStepSize2'
    isvalid = isa(value, 'double') && (length(value) == 1 || length(value) == n);  
  case 'UseVectorized'
    isvalid = isa(value, 'logical');
  case 'LargeScale'
    isvalid = isa(value, 'logical');
  case 'FunValCheck'
    isvalid = isa(value, 'logical'); 
  case 'ValidateInput'
    isvalid = isa(value, 'logical');
  otherwise
    isvalid = true; 
end
end
