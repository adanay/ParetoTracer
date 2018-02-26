% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [c, r, inserted, index] = recover(tree, depth, x, lb, ub, index)
% Inserts the covering cell of a point into the tree as a leaf.
%
% tree is a binary tree that represents a partition of the space. The
% smallest partitions are called cells. Each leaf of the tree represents an
% existent cell.  
% depth is the number of subdivision iterations (depth of the tree).
% Initially make tree = scells.treenode for an empty tree.
% x is the point to insert. x must be a vector. 
% lb and ub are vectors that represent the box constraints.
% index is the index of x in an external list.
%
% Returns the center and radius of the covering cell of x.
% Also returns inserted = true if the covering cell of x was inserted, and
% inserted = false if the cell was not inserted because it was already in
% the tree.
% If the cell was already in the tree, returns the corresponding index in 
% an external list.

x = x(:)';
n = length(x);
inserted = false;
node = tree; % current node

for i = 1 : depth
  j = mod(i - 1, n) + 1; % current subdivision coordinate
    
  % center point at the current subdivision coordinate
  center = lb(j) + (ub(j) - lb(j)) / 2; 
    
  if x(j) < center % the point is at the left side
    if isempty(node.node1)
      node.node1 = scells.treenode; 
      inserted = true; % there was an insertion
    end
        
    node = node.node1;
    ub(j) = center; % deleting the right side       
  else % the point is at the right side
    if isempty(node.node2)
      node.node2 = scells.treenode;
      inserted = true; % there was an insertion
    end
        
    node = node.node2;
    lb(j) = center; % deleting the left side
  end
end

% appends index information
if exist('index', 'var')
  if inserted
    node.index = index;
  end
else
  index = 0;
end

% extracts index information
if ~inserted
  index = node.index;
end

r = (ub - lb) / 2; % radius of the corresponding cell of x
c = lb + r; % center of the corresponding cell of x
end

