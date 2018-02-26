% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Lx1] = bfgsmodchol(Lx0, x0, Jx0, x1, Jx1, opts)
% BFGS update of the Hessians Cholesky factors.
% Bx0 = (Lx0 * Lx0').
% Bx1 = (Lx1 * Lx1').
% Bx0 can be empty. 
%
% Lx0 is the Hessian Cholesky factor to be updated. The size is assumed to  
% be (m x n x n x nobj) where m is the number of invividuals and nobj the  
% no of objectives. In case m = 1 then Lx0 is assumed to have a size 
% of (n x n x nobj).
% x0 and x1 are the previous and current iteration values. For simplicity, 
% if they are vectors, they must be row vectors. Otherwise they will be 
% taken as several individuals of only one variable.
% Jx0 and Jx1 are the previous and current iteration Jacobian values.
% If more than one individual is specified, Jx0 and Jx1 must have a size of 
% (m x nobj x n). If m = 1, they must be of size (nobj x n).
%
% Returns the updated Hessian Cholesky factor Lx1.

[Lx0, x0, Jx0, x1, Jx1, m, n, nobj] = val.valqnmodcholargin(Lx0, x0, Jx0, x1, Jx1, opts);

m0 = size(x0, 1); % m is the number of individuals
m1 = size(x1, 1); % m is the number of individuals

% In case all x0 values are empty, this will be considered the first
% iteration and will return the identity.
if m0 == 0
  if m1 > 1
    Lx1 = utils.eyen([m1 n n nobj], 2, 3); % (m1 x n x n x nobj)
  else
    Lx1 = utils.eyen([n n nobj], 1, 2); % (n x n x nobj)
  end
  return  
end

if m0 == 1
  Jx0 = shiftdim(Jx0, -1); % (1 x nobj x n)
  Lx0 = shiftdim(Lx0, -1); % (1 x n x n x nobj)
  if m1 > m0
    Lx0 = repmat(Lx0, m1, 1, 1, 1);
  end
end
if m1 == 1
  Jx1 = shiftdim(Jx1, -1); % (1 x nobj x n)
end

Rx0 = permute(Lx0, [1 3 2 4]);
Bx0 = eval.lltimes(Lx0);
  
s = bsxfun(@minus, x1, x0); % (m x n) 
y = bsxfun(@minus, Jx1, Jx0); % (m x nobj x n)
    
sB = eval.vhtimes(Bx0, s); % (m x nobj x n) 
sBs = eval.vhvtimes(sB, s); % (m x nobj)
    
sy = eval.vhvtimes(y, s); % (m x nobj)
                 
X = bsxfun(@rdivide, y, sqrt(sy));
X = permute(X, [1 3 4 2]);
R = utils.cholupdaten(Rx0, X, 2, 3);

X = bsxfun(@rdivide, sB, sqrt(sBs));
X = permute(X, [1 3 4 2]);
Rx1 = utils.cholupdaten(R, X, 2, 3, '-');
 
ind = isnan(Rx1) | isinf(Rx1) | ~isreal(Rx1);
Rx1(ind) = Rx0(ind);

Lx1 = permute(Rx1, [1 3 2 4]);

if m == 1
  Lx1 = shiftdim(Lx1, 1); % (n x n x nobj)
end
end

