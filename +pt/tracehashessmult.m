function [hashessmult] = tracehashessmult(multfun)
% Determines whether the PC algorithm is working with Hessian product 
% vectors provided by the user.

% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

hashessmult = ~isempty(multfun) && isstruct(multfun) && (~isempty(multfun.Hwv) || ~isempty(multfun.vH));
end

