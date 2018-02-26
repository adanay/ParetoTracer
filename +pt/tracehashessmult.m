% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [hashessmult] = tracehashessmult(multfun)
% Determines whether the PC algorithm is working with Hessian product 
% vectors provided by the user.

hashessmult = ~isempty(multfun) && isstruct(multfun) && (~isempty(multfun.Hwv) || ~isempty(multfun.vH));
end

