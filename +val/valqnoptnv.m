function [isvalid] = valqnoptnv(name, value)
% Validates the FD options.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if isempty(value)
  isvalid = true;
  return
end

switch(name)
  case 'FunValCheck'
    isvalid = isa(value, 'logical'); 
  case 'ValidateInput'
    isvalid = isa(value, 'logical');
  otherwise
    isvalid = true; 
end
end
