% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [x, v, w, m] = valmopargin(x, v, w, n, nobj, doval)
% Validates the input for a vectorized mop function.

x = x(:, :);
if isvector(x)
  if n == 1
    x = x(:);
  else
    x = x(:)';
  end
end
v = v(:, :);
if isvector(v)
  if n == 1
    v = v(:);
  else
    v = v(:)';
  end
end
w = w(:, :);
if isvector(w)
  if nobj == 1
    w = w(:);
  else
    w = w(:)';
  end
end

if ~doval
  m = max([size(x, 1), size(v, 1), size(w, 1)]);
  return
end

mx = size(x, 1);

if size(x, 2) ~= n || mx < 1
  error('val:valmopargin:InvalidXLength', 'The function is expecting a vector x of length %d.', n);
end

if ~isempty(v)
  if size(v, 2) ~= n
    error('val:valmopargin:InvalidVLength', 'The function is expecting a vector v of length %d.', n);
  end
  mv = size(v, 1);
else
  mv = 0;
end

if ~isempty(w)
  [mw, nobjw] = size(w);
  if nobjw ~= nobj
    error('val:valmopargin:InvalidWLength', 'The function is expecting a vector w of length %d.', nobj);
  end
else
  mw = 0;
end

mm = [mx, mv, mw];
m = max(mm);
mm = mm(mm > 1);
if ~isempty(m) && ~all(mm == m)
  error('val:valmopargin:XAndVWInconsistentM', 'The number of individuals is inconsistent in the specified x and multiply vector.');
end
end

