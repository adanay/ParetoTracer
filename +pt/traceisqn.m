% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [isqn] = traceisqn(objfun, multfun, opts)
% Determines whether the PC algorithm is performing a QN approach.

isqn = ~pt.tracehasfullhess(objfun, multfun) && ~pt.tracehashessmult(multfun) &&...
        strcmpi(opts.HessApprox, 'bfgs');
end

