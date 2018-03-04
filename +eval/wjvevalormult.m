function [wJvx, wJvcount, wJvundef,...
         wJx, wJcount, wJundef,... 
         Jvx, Jvcount, Jvundef,...
         Jx, Jcount, Jundef] =... 
  wjvevalormult(wJv, wJ, Jv, J, x, wJx, Jvx, Jx, w, v, iswarning, forcejac, opts)
% Vectorized Jacobian multiply function evaluation.
% wJv = w' * J * v.
% If wJv is empty, it will perform the multiplication.
% Each row of x is considered an individual to be evaluated.
% Each row of w is considered an individual to be evaluated.
% Each row of v is considered an individual to be evaluated.
% x, w, and v must have the same no of individuals or only one (in order to be repeated).
% iswarning defines whether to display a warning or an error if the
% Jacobian value is undefined. The Jacobian value will be checked only if 
% opts.FunValCheck is true.
% forcejac will force one evaluation of the Jacobian to be returned in Jx 
% even if wJv (or wJ or Jv) was passed and evaluated.
% opts are the optimization options.
%
% The result is always of size (m x 1).
% Jx (if entered) is assumed to be of size (m x obj x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% wJx (if entered) is assumed to be of size (m x n) even if m = 1 
% (after being calculated using one of the vec eval functions).
% Jvx (if entered) is assumed to be of size (m x nobj) even if m = 1 
% (after being calculated using one of the vec eval functions).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

forcejac = forcejac && nargout > 9;

wJvcount = 0;
wJcount = 0;
Jvcount = 0;
Jcount = 0;

wJvundef = false;
wJundef = false;
Jvundef = false;
Jundef = false;

n = size(x, 2); % n is the number of variables
nobj = size(w, 2); % nobj is the number of objectives

if isempty(wJv)
  if isempty(wJ) && isempty(wJx) && isempty(Jv) && isempty(Jvx)
    wJv_J();
  else
    if nobj > n
      if isempty(wJx) % favor wJ
        if isempty(wJ)
          wJv_Jv();
        else
          if isempty(Jvx)
            wJv_wJ();
          else
            wJv_Jv();
          end
        end
      else
        wJv_wJ();
      end
    else % favor Jv 
      if isempty(Jvx)
        if isempty(Jv)
          wJv_wJ();
        else
          if isempty(wJx)
            wJv_Jv();
          else
            wJv_wJ();
          end
        end
      else
        wJv_Jv();
      end
    end
  end
else
  [wJvx, wJvcount, wJvundef] = eval.wjveval(wJv, x, w, v, iswarning, opts);
end

if forcejac && isempty(Jx)
  [Jx, Jcount, Jundef] = eval.jeval(J, x, iswarning, opts);
end

  function [] = wJv_J()
    if isempty(Jx)
      [Jx, Jcount, Jundef] = eval.jeval(J, x, iswarning, opts);
    end
      [wJx, wJcount] = eval.wjtimes(Jx, w);  
      [wJvx, wJvcount] = eval.wjvtimes(wJx, v);
  end

  function [] = wJv_wJ()
    if isempty(wJx)
      [wJx, wJcount, wJundef] = eval.wjeval(wJ, x, w, iswarning, opts); % (m x n)
    end
    [wJvx, wJvcount] = eval.wjvtimes(wJx, v);
  end

  function [] = wJv_Jv()
    if isempty(Jvx)
      [Jvx, Jvcount, Jvundef] = eval.jveval(Jv, x, v, iswarning, opts); % (m x nobj)
    end
    [wJvx, wJvcount] = eval.jvwtimes(Jvx, w);
  end
end