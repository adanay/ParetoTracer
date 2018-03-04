function [issd] = minissd(objfun, multfun, opts)
% Determines whether the optimization algorithm is completely excluding the 
% Hessian information, i.e., an steepest descent approach.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

ishessprov = isstruct(objfun) && ~isempty(objfun.H) ||...
             ~isempty(multfun) && isstruct(multfun) && ~isempty(multfun.vH); 
        
issd = ~ishessprov && strcmpi(opts.HessApprox, 'off');
end

