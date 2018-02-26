% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [Bx1] = bfgs(Bx0, x0, Jx0, x1, Jx1, opts)
% BFGS update of the Hessians.
%
% Bx0 is the Hessian to be updated. The size is assumed to be 
% (m x n x n x nobj) where m is the number of invividuals and nobj the no 
% of objectives. In case m = 1 then Bx0 is assumed to have a size 
% of (n x n x nobj).
% x0 and x1 are the previous and current iteration values. For simplicity, 
% if they are vectors, they must be row vectors. Otherwise they will be 
% taken as several individuals of only one variable.
% Jx0 and Jx1 are the previous and current iteration Jacobian values.
% If more than one individual is specified, Jx0 and Jx1 must have a size of 
% (m x nobj x n). If m = 1, they must be of size (nobj x n).
%
% Returns the updated Hessian Bx1.

[Bx0, x0, Jx0, x1, Jx1, m, n, nobj] = val.valqnargin(Bx0, x0, Jx0, x1, Jx1, opts);

m0 = size(x0, 1); % m is the number of individuals
m1 = size(x1, 1); % m is the number of individuals

% In case all x0 values are empty, this will be considered the first
% iteration and will return the identity.
if m0 == 0
  if m1 > 1
    Bx1 = utils.eyen([m1 n n nobj], 2, 3); % (m1 x n x n x nobj)
  else
    Bx1 = utils.eyen([n n nobj], 1, 2); % (n x n x nobj)
  end
  return 
end

if m0 == 1
  Jx0 = shiftdim(Jx0, -1); % (1 x nobj x n)
  Bx0 = shiftdim(Bx0, -1); % (1 x n x n x nobj)
  if m1 > m0
    Bx0 = repmat(Bx0, m1, 1, 1, 1);
  end
end
if m1 == 1
  Jx1 = shiftdim(Jx1, -1); % (1 x nobj x n)
end
  
s = bsxfun(@minus, x1, x0); % (m x n)
y = bsxfun(@minus, Jx1, Jx0); % (m x nobj x n)
  
sB = eval.vhtimes(Bx0, s); % (m x nobj x n)  

BssB = eval.jjtimes(sB); % (m x n x n x nobj)

sBs = eval.vhvtimes(sB, s); % (m x nobj)
sBs = permute(sBs, [1 3 4 2]); % (m x 1 x 1 nobj);
    
yy = eval.jjtimes(y); % (m x n x n x nobj)

sy = eval.vhvtimes(y, s); % (m x nobj)
sy = permute(sy, [1 3 4 2]); % (m x 1 x 1 nobj);
  
Bx1 = bsxfun(@minus, Bx0, bsxfun(@rdivide, BssB, sBs)) +...
                          bsxfun(@rdivide, yy, sy); % BFGS update
 
ind = isnan(Bx1) | isinf(Bx1) | ~isreal(Bx1);
Bx1(ind) = Bx0(ind);

if m == 1
  Bx1 = shiftdim(Bx1, 1); % (n x n x nobj)
end
end
