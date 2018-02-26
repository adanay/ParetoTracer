% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [vHx, vHcount, vHundef, vHapprox,... 
          Hx, Hcount, Hundef, Happrox, Hident,... 
          fx, fcount, Jx, Jcount, Jundef] =...
  vhevalormultorfd(vH, H, f, J, x, fx, Jx, Hx, v, lb, ub, iswarning, forcehess, formident, opts)
% Vectorized Hessian multiply function evaluation.
% vH = [v' * H1; v' * H2; ...; v' * Hnobj].
% If vH is empty, it will perform the multiplication.
% If H is also empty, it will use finite differences (FD).
% Each row of x is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and vH(x) or 
% H(x) is undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be computed by FD if the Jacobian or the objective functions 
% are available. Otherwise it will be approximated to the identity.
% + iswarning = false implies that an error will be thrown.
% In case that FD is applied with the Jacobian, the value of the Jacobian 
% is the one that will be checked if opts.FunValCheck is true. If also
% iswarning = true, then FD will be re-computed using the objective 
% function f if available.
% forcehess will force one evaluation of the Hessian to be returned in Hx 
% even if vH was passed and evaluated.
% formident will force the formation of the identity to be returned in Hx 
% in case there is no other choice than approximating Hx to the identity. 
% Otherwise Hx will be empty (but known to be the identity if Hident is true).
% opts are the optimization options.
%
% The result is always of size (m x nobj x n) even if m = 1.
% fx (if entered) is assumed to be of size (m x obj) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Hx (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).

forcehess = forcehess && nargout > 4;

vHcount = 0;
vHundef = false;
vHapprox = false;

Hcount = 0;
Hundef = false;
Happrox = false;
Hident = false;

fcount = 0;

Jcount = 0;
Jundef = false;

mx = size(x, 1); % m is the number of individuals
mv = size(v, 1); 
m = max(mx, mv); 

if isempty(vH) && isempty(H) && isempty(Hx)
  vhess();
else
  [vHx, vHcount, vHundef, Hx, Hcount, Hundef] = eval.vhevalormult(vH, H, x, Hx, v, iswarning, forcehess, opts);
  if vHundef || Hundef
    if ~isempty(f) || ~isempty(J)
      Hx = [];
      vhess();
    else
      nobj = size(vHx, 2);
      vHapprox = true;
      Happrox = true;
      Hident = true;
      [vHx, ~, ~, ~, Hx] = eval.vhevalormultorsd([], [], x, [], v, nobj, iswarning, forcehess, formident, opts);
    end
  end
end

  function [] = vhess()
    vHapprox = true;
    
    if forcehess
      vhess_H();
    else
      if mx == 1 && ~isempty(Jx)
        Jx = shiftdim(Jx, 1);
      end
      
      try
        [vHx, fx, fcount, Jx, Jcount, Jundef] = fd.hessian(f, J, x, fx, Jx, lb, ub, v, opts);
        if m == 1
          vHx = shiftdim(vHx, -1); % (1 x nobj x n)
        end
      catch % When a fixed direction v is given for FD, 
            % it is possible that x + c * v is not feasible no matter how small c > 0 is. 
        vhess_H();
      end    
      
      if Jundef
        if ~isempty(J) && ~isempty(f)
          try
            [vHx, fx, fcount] = fd.hessian(f, [], x, fx, [], lb, ub, v, opts);
            if m == 1
              vHx = shiftdim(vHx, -1); % (1 x nobj x n)
            end
          catch % When a fixed direction v is given for FD, 
                % it is possible that x + c * v is not feasible no matter how small c > 0 is.
            vhess_H();
          end
        else
          nobj = size(vHx, 2);
          Happrox = true;
          Hident = true;
          vHx = resizev(v, nobj);
        end
      end
      
      if mx == 1 && ~isempty(Jx)
        Jx = shiftdim(Jx, -1); % (1 x nobj x n)
      end
    end
  end

  function [] = vhess_H()
    if isempty(Hx)
      [Hx, ~, ~, Happrox, Hident, fx, fcount, Jx, Jcount, Jundef, nobj] = eval.hevalorfd([], f, J, x, fx, Jx, lb, ub, iswarning, forcehess && formident, opts);
    end
    if Hident
      vHx = resizev(v, nobj);
    else
      [vHx, tvHcount] = eval.vhevalormult([], [], x, Hx, v, iswarning, false, opts);
      vHcount  = vHcount + tvHcount;
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