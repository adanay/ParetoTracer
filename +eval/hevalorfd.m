function [Hx, Hcount, Hundef, Happrox, Hident,...
          fx, fcount, Jx, Jcount, Jundef, nobj] = hevalorfd(H, f, J, x, fx, Jx, lb, ub, iswarning, formident, opts)
% Vectorized Hessian function evaluation.
% If H is empty, finite differences (FD) are used.
% Each row of x is considered an individual to be evaluated.
% iswarning defines whether to display a warning or an error if the
% Hessian value is undefined. The Hessian value will be checked only if 
% opts.FunValCheck is true. I.e., if opts.FunValCheck is true and H(x) is 
% undef, then
% + iswarning = true implies that a warning will be displayed and the
% Hessian will be computed by FD if the Jacobian or the objective functions 
% are available. Otherwise the Hessian will be approximated to the identity 
% matrix.
% + iswarning = false implies that an error will be thrown.
% In case that FD is applied with the Jacobian, the value of the Jacobian 
% is the one that will be checked if opts.FunValCheck is true. If also
% iswarning = true, then FD will be re-computed using the objective 
% function f if available.
% formident will force the formation of the identity to be returned in Hx 
% in case there is no other choice than approximating Hx to the identity. 
% Otherwise Hx will be empty (but known to be the identity if Hident is true).
% opts are the optimization options.
%
% The result is always of size (m x n x n x nobj) even if m = 1.
% The result of H(x) is assumed to be (n x n x nobj) if m = 1.
% fx (if entered) is assumed to be of size (m x obj) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

Hcount = 0;
Hundef = false;
Happrox = false;
Hident = false; % the Hessian was approx to the identity matrix

fcount = 0;

Jcount = 0;
Jundef = false;

m = size(x, 1); % m is the number of individuals

if isempty(H) 
  hess();
else
  [Hx, Hcount, Hundef] = eval.heval(H, x, iswarning, opts);
  nobj = size(Hx, 4);
  if Hundef
    if ~isempty(f) || ~isempty(J)
      hess();
    else
      Happrox = true;
      Hident = true;
      Hx = eval.hevalorsd([], x, nobj, iswarning, formident, opts);
    end
  end
end

  function [] = hess()
    Happrox = true;
    
    if m == 1 && ~isempty(Jx)
      Jx = shiftdim(Jx, 1);
    end
    
    [Hx, fx, fcount, Jx, Jcount, Jundef] = fd.hessian(f, J, x, fx, Jx, lb, ub, [], opts);
    nobj = size(Hx, 4);
    if m == 1
      Hx = shiftdim(Hx, -1); % (1 x n x n x nobj)
    end
    if Jundef
      if ~isempty(f)
        [Hx, fx, fcount] = fd.hessian(f, [], x, fx, [], lb, ub, [], opts);
        if m == 1
          Hx = shiftdim(Hx, -1); % (1 x n x n x nobj)
        end
      else
        Hident = true;
        Hx = eval.hevalorsd([], x, nobj, iswarning, formident, opts);
      end
    end
    
    if m == 1 && ~isempty(Jx)
      Jx = shiftdim(Jx, -1); % (1 x nobj x n)
    end    
  end
end