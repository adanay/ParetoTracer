function [Hwvx, Hwvcount, Hwvundef,... 
          Hwx, Hwcount, Hwundef,... 
          vHx, vHcount, vHundef,... 
          Hx, Hcount, Hundef] =...
  hwvevalormult(Hwv, Hw, vH, H, x, Hwx, vHx, Hx, w, v, iswarning, forcehess, opts)
% Vectorized Hessian multiply function evaluation.
% Hwv = (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v.
% If Hwv is empty, it will perform the multiplication.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x, w, and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true.
% forcehess will force one evaluation of the Hessian to be returned in Hx 
% even if Hwv (or Hw or vH) was passed and evaluated.
% opts are the optimization options.
%
% The result is always of size (m x n) even if m = 1.
% Hx (if entered) is assumed to be of size (m x n x n x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Hwx (if entered) is assumed to be of size (m x n x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% vHx (if entered) is assumed to be of size (m x nobj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcehess = forcehess && nargout > 9;

Hwvcount = 0;
Hwvundef = false;

Hwcount = 0;
Hwundef = false;

vHcount = 0;
vHundef = false;

Hcount = 0;
Hundef = false;

n = size(x, 2); % n is the number of variables
nobj = size(w, 2); % nobj is the number of objectives

if isempty(Hwv)
  if isempty(Hw) && isempty(Hwx) && isempty(vH) && isempty(vHx)
    Hwv_H();
  else
    if nobj > n % favor Hw (n x n)
      if isempty(Hwx)
        if isempty(Hw)
          Hwv_vH(); % better vH (nobj x n) than H (n x n x nobj)
        else
          if isempty(vHx)
            Hwv_Hw();
          else
            Hwv_vH(); % vH (nobj x n) already computed so no need to evaluate Hw
          end
        end
      else
        Hwv_Hw();
      end
    else % favor vH
      if isempty(vHx)
        if isempty(vH)
          Hwv_Hw(); % better Hw (n x n) than H (n x n x nobj)
        else
          if isempty(Hwx)
            Hwv_vH();
          else
            Hwv_Hw(); % Hw (n x n) already computed so no need to evaluate vH
          end
        end
      else
        Hwv_vH();
      end
    end
  end
else
  [Hwvx, Hwvcount, Hwvundef] = eval.hwveval(Hwv, x, w, v, iswarning, opts);
end

if forcehess && isempty(Hx)
  [Hx, Hcount, Hundef] = eval.heval(H, x, iswarning, opts);
end

  function [] = Hwv_H()
    if isempty(Hx)
      [Hx, Hcount, Hundef] = eval.heval(H, x, iswarning, opts);
    end
    [Hwx, Hwcount] = eval.hwtimes(Hx, w);
    [Hwvx, Hwvcount] = eval.hwvtimes(Hwx, v);
  end

  function [] = Hwv_Hw()
    if isempty(Hwx)
      [Hwx, Hwcount, Hwundef] = eval.hweval(Hw, x, w, iswarning, opts); % (m x n x n)
    end
    [Hwvx, Hwvcount] = eval.hwvtimes(Hwx, v);
  end

  function [] = Hwv_vH()
    if isempty(vHx)
      [vHx, vHcount, vHundef] = eval.vheval(vH, x, v, iswarning, opts); % (m x nobj x n)
    end
    [Hwvx, Hwvcount] = eval.vhwtimes(vHx, w);
  end
end