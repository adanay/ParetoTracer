function [hasfullhess] = tracehasfullhess(objfun, multfun)
% Determines whether the PC algorithm is working with full Hessians 
% provided by the user.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

hasfullhess = isstruct(objfun) && ~isempty(objfun.H) ||...
             ~isempty(multfun) && isstruct(multfun) && ~isempty(multfun.Hw);
end

