% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [opts] = defopts()
% Default options for Quasi-Newton (QN) updates. 

opts = struct(...
  'FunValCheck', true,...
  'ValidateInput', true);
end
