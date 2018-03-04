function [Hx1, Hcount, Hundef, Happrox, Hident, Lx1, Lmodif] = hevalorqnandmodchol(H, x1, Jx1, x0, Jx0, Lx0, iswarning, formident, opts)
% Vectorized Hessian function evaluation plus a modified Cholesky 
% decomposition.
% If H is empty, a modified Cholesky QN update will be performed and stored 
% in Lx1. In this case, the returned Hx1 is empty. (It can be later 
% obtained as Hx1 = Lx1 * Lx1'.) 
% If, on the other hand, H is provided, the Hessian will be computed and  
% stored in Hx1. The Cholesky modification of Hx1 will be stored in Lx1. 
% If Lmodif = true, then Hx1 is not positive definite, thus Lx1 * Lx1'   
% needs to be computed.
% Each row of x1, (and x0) is considered an individual to be evaluated.
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and H(x1) is 
% undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be updated if the previous iteration values are provided.
% Otherwise it will be approximated to the identity.
% + iswarning = false implies that an error will be thrown.
% formident will force the formation of the identity to be returned in Hx1
% (and/or Lx1) in case there is no other choice than approximating Hx1 to  
% the identity. This may happen e.g. if the previous iteration values are  
% not provided. Otherwise Hx1 (and Lx1) will be empty (but known to be the 
% identity if Hident is true).
% opts are the optimization options.
%
% The result is always of size (m x n x n x nobj) even if m = 1.
% The result of H(x) is assumed to be (n x n x nobj) if m = 1.
% Jx1 (and Jx0) (if entered) is assumed to be of size (m x obj x n) even  
% if m = 1 (after being calculated using one of the vec eval functions).
% Hx0 (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

Hx1 = [];
Hcount = 0;
Hundef = false;
Happrox = false;
Hident = false;

Lx1 = [];
Lmodif = [];

m0 = size(x0, 1); 
m1 = size(x1, 1); 
m = max(m0, m1); % m is the number of individuals

if isempty(H)                
  hess();
else
  [Hx1, Hcount, Hundef] = eval.heval(H, x1, iswarning, opts);
  if Hundef
    if ~isempty(Jx1) && ~isempty(x0) && ~isempty(Jx0)
      Hx1 = [];
      hess();
    else
      nobj = size(Hx1, 4);
      Happrox = true;
      Hident = true;
      
      Hx1 = eval.hevalorsd([], x1, nobj, iswarning, formident, opts);
      
      if nargout > 5
        Lmodif = false;
        if formident
          Lx1 = Hx1;
        end
      end
    end
  else
    if nargout > 5
      [Lx1, Lmodif] = utils.modcholn(Hx1, 2, 3);
    end
  end
end

  function [] = hess()
    Happrox = true;
    Lmodif = true;
    
    if m0 == 1
      Jx0 = shiftdim(Jx0, 1);
      Lx0 = shiftdim(Lx0, 1);
    end
    if m1 == 1
      Jx1 = shiftdim(Jx1, 1);
    end
    
    Lx1 = qn.bfgsmodchol(Lx0, x0, Jx0, x1, Jx1, opts);
    
    if m == 1
      Lx1 = shiftdim(Lx1, -1); % (1 x n x n x nobj)
    end
  end
end