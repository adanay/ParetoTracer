function [c, r, index] = ballcontains(tree, depth, list, x, e, lb, ub)
% Attempts to determine whether there already exists a point close to x. It  
% is important to note that only points contained in the neighboring cells  
% of x will be considered, i.e., it is assumed e < r, where r is the radius   
% of the cells. If a neighbor point of x belongs to the e-ball of x, its  
% corresponding cell will be returned (if there is a point in the same cell 
% of x, the cell of x will be returned as in scells.treecontains). 
% Note that since the number of neighbor cells grow exponentially with n,
% the algorithm may consider to check only neighbors in coordinate 
% directions.
%
% tree is a binary tree that represents a partition of the space. The
% smallest partitions are called cells. Each leaf of the tree represents an
% existent cell.
% depth is the number of subdivision iterations (depth of the tree).
% list is an external list that contains the existent points (whose 
% covering cells are represented by a leaf of the tree).
% x is the point to consider whether or not is alone in the neighborhood.
% e is the neighborhood radius.
% lb and ub are vectors that represent the box constraints.
%
% Returns the center and radius of the covering cell of x (or a close 
% neighbor of x).

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

x = x(:)';

% check first if there is a point in the same cell of x (conventional)
[c, r] = scells.treecontains(tree, depth, x, lb, ub);
if ~isempty(c)
  return
end

% covering cell of x
[c, r] = scells.cellof(x, lb, ub, depth);
if isempty(c)
  return
end

% bounds of the covering cell of x
l = c - r;
u = c + r;

% active bounds: those that intercepts the neighborhood
la = find(x - l < e);
r1 = length(la);

ua = find(u - x < e);
r2 = length(ua);

a = [-la, ua];

if r1 + r2 > 6 % O((r1 + r2))
  for i = 1 : r1 + r2 
    % verifies whether there is a close neighbor of x in the neighbor 
    % cell defined by the chosen dimensions
    [cn, rn, index] = ballcontains(a(i));
    if ~isempty(cn)
      c = cn;
      r = rn;
      return
    end
  end
else % O(2^(r1 + r2))
  for i = 1 : r1 + r2 
    dims = nchoosek(a, i); % all possible combinations of i dimensions 

    for j = 1 : size(dims, 1) % number of neighbor cells to be analized
      dim = dims(j, :);

      % check for repeated dimensions, i.e., it is imposible to find a cell
      % that is a neighbor of this set of dimensions simultaneously.
      if length(unique(abs(dim))) ~= length(dim)
        continue
      end

      % verifies whether there is a close neighbor of x in the neighbor 
      % cell defined by the chosen dimensions
      [cn, rn, index] = ballcontains(dim);
      if ~isempty(cn)
        c = cn;
        r = rn;
        return
      end
    end
  end
end

% nothing was found
c = [];
r = [];
index = 0;

function [cn, rn, index] = ballcontains(dim)
  % center of the neighbor cell
  cn = c;
  cn(abs(dim)) = c(abs(dim)) + sign(dim) .* 2 .* r(abs(dim));

  [cn, rn, index] = scells.treecontains(tree, depth, cn, lb, ub);

  if ~isempty(cn) && index ~= 0 && norm(x - list(index, :)) <= e
    return
  else
    cn = [];
    rn = [];
    index = 0;
  end
end
end

