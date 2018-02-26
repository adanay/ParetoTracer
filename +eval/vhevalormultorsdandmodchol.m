% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [vHx, vHcount, vHundef, vHapprox,...
          Hx, Hcount, Hundef, Happrox, Lx, Lmodif] = vhevalormultorsdandmodchol(vH, H, x, Hx, v, nobj, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation plus a modified Cholesky
% decomposition.
% vH = [v' * H1; v' * H2; ...; v' * Hnobj], where H1,..., Hnobj are modified
% to ensure they are positive definite.
% If vH is empty, it will perform the multiplication.
% If H is also empty, the Hessian will be approximated to the identity matrix.
% If Lmodif = true, then Hx is not positive definite, thus Lx * Lx' needs  
% to be computed.
% If Lmodif is empty, then the modification did not take place as e.g. the
% multiply function vH was provided or Hx was provided.
% Each row of x is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and vH(x) or 
% H(x) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be approximated to the identity matrix.
% + iswarning = false implies that an error will be thrown.
% forcehess will force one evaluation of the Hessian to be returned in Hx 
% even if vH was passed and evaluated.
% formident will force the formation of the identity to be returned in Hx 
% (and Lx) if it is approximated. Otherwise Hx (and Lx) will be empty (but 
% known to be the identity if Happrox is true).
% The Cholesky modification will be performed only if the Hessian is
% defined and if it had to be computed. If the Hessian is passed as an
% argument it will be assumed to be positive definite and will not be
% modified.
% opts are the optimization options.
%
% The result is always of size (m x nobj x n) even if m = 1.
% Hx (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).

forcehess = forcehess && nargout > 4;

vHcount = 0;
vHundef = false;
vHapprox = false;

Hcount = 0;
Hundef = false;
Happrox = false;

Lx = [];
Lmodif = [];

[mx, n] = size(x); % m is the number of individuals and n is the number of variables.
mv = size(v, 1); 
m = max(mx, mv);

if isempty(vH) && isempty(H) && isempty(Hx)
  vhess(nargout);
else
  if nargout > 8
    [vHx, vHcount, vHundef, Hx, Hcount, Hundef, Lx, Lmodif] = eval.vhevalormultandmodchol(vH, H, x, Hx, v, iswarning, forcehess, opts);
  else
    [vHx, vHcount, vHundef, Hx, Hcount, Hundef] = eval.vhevalormultandmodchol(vH, H, x, Hx, v, iswarning, forcehess, opts);
  end
  if vHundef || Hundef
    Hx = [];
    vhess(nargout);
  end
end

  function [] = vhess(onargout)
    vHapprox = true;
    Happrox = true;
    
    vHx = resizev(v, nobj);
    
    if forcehess && isempty(Hx)  
      if formident
        Hx = utils.eyen([m n n nobj], 2, 3);
      end
      
      if onargout > 8
        Lmodif = false;
        if formident
          Lx = Hx;
        end
      end
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