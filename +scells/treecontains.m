% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [c, r, index] = treecontains(tree, depth, x, lb, ub)
% Determines whether the covering cell of a point is into the tree.
% 
% tree is a binary tree that represents a partition of the space. The
% smallest partitions are called cells. Each leaf of the tree represents an
% existent cell.  
% depth is the number of subdivision iterations (depth of the tree).
% x is the point to verify. x must be a vector. 
% lb and ub are vectors that represent the box constraints.
%
% Returns the center and radius of the covering cell of x if it exists.
% If contained, also returns the associated index of x in an external list 
% if exists.

c = [];
r = [];
index = 0;

if ~utils.isfeasible(x, lb, ub)
  return
end

n = length(x); % number of variables
node = tree; % current node

for i = 1 : depth
  j = mod(i - 1, n) + 1; % current subdivision coordinate
    
  % center point at the current subdivision coordinate
  center = lb(j) + (ub(j) - lb(j)) / 2; 
    
  if x(j) < center % the point is at the left side
    if isempty(node.node1)
      return
    end
        
    node = node.node1;
    ub(j) = center; % deleting the right side       
  else % the point is at the right side
    if isempty(node.node2)
      return
    end
        
    node = node.node2;
    lb(j) = center; % deleting the left side
  end
end

index = node.index;
if isempty(index)
  index = 0;
end

r = (ub - lb) / 2; % radius of the corresponding cell of x
c = lb + r; % center of the corresponding cell of x
end