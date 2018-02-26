% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [undef] = checkval(x, name, iswarning)
% Determines whether a value has an undefined entry an throws an  
% error if it is the case.
% If iswarning = true, a warning will be displayed instead.

undef = false;

if nargin < 3
  iswarning = false;
end

if ~isa(x, 'double')
  undef = true;
  if iswarning
    warning('val:checkval:NonDouble', 'The value of %s is not numeric.', name);
  else
    error('val:checkval:NonDouble', 'The value of %s is not numeric.', name);
  end
end

if ~isreal(x)
  undef = true;
  if iswarning
    warning('val:checkval:Complex', 'The value of %s has a nonreal entry.', name);
  else
    error('val:checkval:Complex', 'The value of %s has a nonreal entry.', name);
  end
end

if any(isinf(x(:)))
  undef = true;
  if iswarning
    warning('val:checkval:Inf', 'The value of %s has an Inf entry.', name);
  else
    error('val:checkval:Inf', 'The value of %s has an Inf entry.', name);
  end
end

if any(isnan(x(:)))
  undef = true;
  if iswarning
    warning('val:checkval:NaN', 'The value of %s has a NaN entry.', name);
  else
    error('val:checkval:NaN', 'The value of %s has a NaN entry.', name);
  end
end
end
