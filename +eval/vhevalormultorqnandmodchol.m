function [vHx1, vHcount, vHundef, vHapprox,...
          Hx1, Hcount, Hundef, Happrox, Hident, Lx1, Lmodif] = vhevalormultorqnandmodchol(vH, H, x1, Jx1, Hx1, x0, Jx0, Lx0, v, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation plus a modified Cholesky
% decomposition.
% vH = [v' * H1; v' * H2; ...; v' * Hnobj], where H1,..., Hnobj are modified
% to ensure they are positive definite.
% If vH is empty, it will perform the multiplication.
% If H is also empty, a modified Cholesky QN update will be performed and  
% stored in Lx1. In this case, the returned Hx1 is empty. (It can be later 
% obtained as Hx1 = Lx1 * Lx1'.) 
% If Lmodif = true, then Hx is not positive definite, thus Lx * Lx' needs  
% to be computed.
% If Lmodif is empty, then the modification did not take place as e.g. the
% multiply function vH was provided or Hx was provided.
% Each row of x1 (and x0) is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x1 and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and vH(x1) or 
% H(x1) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be updated if the previous iteration values are provided.
% Otherwise it will be approximated to the identity matrix.
% + iswarning = false implies that an error will be thrown.
% forcehess will force one evaluation of the Hessian to be returned in Hx1 
% even if vH was passed and evaluated.
% formident will force the formation of the identity to be returned in Hx1
% (and/or Lx1) in case there is no other choice than approximating Hx1 to  
% the identity. This may happen e.g. if the previous iteration values are  
% not provided. Otherwise Hx1 (and Lx1) will be empty (but known to be the 
% identity if Hident is true).
% The Cholesky modification will be performed only if the Hessian is
% defined and if it had to be computed. If the Hessian is passed as an
% argument it will be assumed to be positive definite and will not be
% modified.
% opts are the optimization options.
%
% The result is always of size (m x nobj x n) even if m = 1.
% Jx1 (and Jx0) (if entered) is assumed to be of size (m x obj x n) even  
% if m = 1 (after being calculated using one of the vec eval functions).
% Hx1 (and Hx0) (if entered) is assumed to be of size (m x n x n x nobj)  
% even if m = 1 (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcehess = forcehess && nargout > 4;

vHcount = 0;
vHundef = false;
vHapprox = false;

Hcount = 0;
Hundef = false;
Happrox = false;
Hident = false;

Lx1 = [];
Lmodif = [];

[mx1, n] = size(x1); % m is the number of individuals and n is the number of variables
mx0 = size(x0, 1); % m is the number of individuals
mx = max(mx0, mx1); % m is the number of individuals
mv = size(v, 1); % m is the number of individuals
m = max(mx, mv); % m is the number of individuals

if isempty(vH) && isempty(H) && isempty(Hx1)
  vhess();
else
  if nargout > 9
    [vHx1, vHcount, vHundef, Hx1, Hcount, Hundef, Lx1, Lmodif] = eval.vhevalormultandmodchol(vH, H, x1, Hx1, v, iswarning, forcehess, opts);
  else
    [vHx1, vHcount, vHundef, Hx1, Hcount, Hundef] = eval.vhevalormultandmodchol(vH, H, x1, Hx1, v, iswarning, forcehess, opts);
  end
  if vHundef || Hundef
    if ~isempty(Jx1) && ~isempty(x0) && ~isempty(Jx0)
      vhess();
    else
      nobj = size(vHx1, 2);
      vHapprox = true;
      Happrox = true;
      Hident = true;
      if nargout > 9
        [vHx1, ~, ~, ~, Hx1, ~, ~, ~, Lx1, Lmodif] = eval.vhevalormultorsdandmodchol([], [], x1, [], v, nobj, iswarning, forcehess, formident, opts);
      else
        [vHx1, ~, ~, ~, Hx1] = eval.vhevalormultorsdandmodchol([], [], x1, [], v, nobj, iswarning, forcehess, formident, opts);
      end
    end
  end
end

  function [] = vhess()
    vHapprox = true;
    Happrox = true;
    
    [Hx1, ~, ~, ~, Hident, Lx1, Lmodif] = eval.hevalorqnandmodchol([], x1, Jx1, x0, Jx0, Lx0, iswarning, forcehess && formident, opts);
    if Hident
      nobj = size(Jx1, 2);
      vHx1 = resizev(v, nobj);
    else
      [vHx1, vHcount] = eval.vlltimes(Lx1, v);
    end
  end

  function[vHx] = resizev(v, nobj)
    vHx = v; % (mv x n)
    vHx = permute(vHx, [1 3 2]); % (mv x 1 x n)
    vHx = repmat(vHx, 1, nobj, 1); % (mv x nobj x n)
    if mv == 1 && mx ~= 1
      vHx = repmat(vHx, m, 1, 1); % (m x nobj x n)
    end
  end
end