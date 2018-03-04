function [vHx1, vHcount, vHundef, vHapprox,...
          Hx1, Hcount, Hundef, Happrox, Hident] = vhevalormultorqn(vH, H, x1, Jx1, Hx1, x0, Jx0, Hx0, v, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation.
% vH = [v' * H1; v' * H2; ...; v' * Hnobj].
% If vH is empty, it will perform the multiplication.
% If H is also empty, a QN update will be performed.
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
% in case there is no other choice than approximating Hx1 to the identity. 
% This may happen e.g. if the previous iteration values are not provided. 
% Otherwise Hx1 will be empty (but known to be the identity if Hident is true).
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

mx = size(x1, 1); % m is the number of individuals
mv = size(v, 1); 
m = max(mx, mv); 

if isempty(vH) && isempty(H) && isempty(Hx1)
  vhess();
else
  [vHx1, vHcount, vHundef, Hx1, Hcount, Hundef] = eval.vhevalormult(vH, H, x1, Hx1, v, iswarning, forcehess, opts);
  if vHundef || Hundef
    if ~isempty(Jx1) && ~isempty(x0) && ~isempty(Jx0)
      vhess();
    else
      nobj = size(vHx1, 2);
      vHapprox = true;
      Happrox = true;
      Hident = true;
      [vHx1, ~, ~, ~, Hx1] = eval.vhevalormultorsd([], [], x1, [], v, nobj, iswarning, forcehess, formident, opts);
    end
  end
end

  function [] = vhess()
    vHapprox = true;
    Happrox = true;
    
    [Hx1, ~, ~, ~, Hident] = eval.hevalorqn([], x1, Jx1, x0, Jx0, Hx0, iswarning, forcehess && formident, opts);
    if Hident
      nobj = size(Jx1, 2);
      vHx1 = resizev(v, nobj);
    else
      [vHx1, tvHcount] = eval.vhevalormult([], [], x1, Hx1, v, iswarning, false, opts);
      vHcount = vHcount + tvHcount;
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