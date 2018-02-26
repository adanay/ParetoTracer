% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [isvalid] = valqnoptnv(name, value)
% Validates the FD options.

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
