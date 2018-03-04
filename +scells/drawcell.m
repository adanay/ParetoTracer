function [h] = drawcell(c, r, varargin)
% Draws a cell with the specified center and radius.
%
% c is the center of the cell.
% r is the radius of the cell.
% Use property/value pairs to specify additional patch properties.
%
% Returns the handle to the graphic object.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if isempty(c) || isempty(r) 
  return
end

n = length(c); % number of variables

x = ([0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0] - 0.5) * 2 * r(1) + c(1);
y = ([0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1] - 0.5) * 2 * r(2) + c(2);

if n == 2
  h = patch(x, y, 'w', 'facealpha', 0.5);
else
  view(3);
  z = ([0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1] - 0.5) * 2 * r(3) + c(3);
  h = patch(x, y, z, 'w', 'facealpha', 0.5);
end

if ~isempty(varargin)
  set(h, varargin{:});
end
end