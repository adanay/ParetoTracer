function [issd] = traceissd(objfun, multfun, opts)
% Determines whether the PC algorithm is completely excluding the 
% Hessian information, i.e., an steepest descent approach.

% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

issd = ~pt.tracehasfullhess(objfun, multfun) && ~pt.tracehashessmult(multfun) &&...
        strcmpi(opts.HessApprox, 'off');
end

