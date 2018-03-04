function [] = printtree(tree, lb, ub)
% Prints the tree built by the subdivision algorithm.
%
% tree is a binary tree that represents a partition of the space 
% constrained by lb and ub. The smallest partitions are called cells. Each 
% leaf of the tree represents an existent cell.  

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

fprintf('\n');
printnode(tree, lb, ub, 0);
end

function [] = printnode(node, lb, ub, depth)
% Prints a node of the tree built by the subdivision algorithm.
% lb and ub are vectors that represent the box constraints.
% depth is the current iteration of the subdivision (the depth of the tree).

if isempty(node)
  return
end

% printing current node
if depth > 0
  fprintf('%d', depth);
end
for i = 1 : depth
  fprintf('-');
end
utils.printdomain(lb, ub);
fprintf('\n');

n = length(lb); % number of variables
j = mod(depth, n) + 1; % current subdivision coordinate

% set subdivision
lb1 = lb;
ub1 = ub;
ub1(j) = lb1(j) + (ub1(j) - lb1(j)) / 2; 

lb2 = lb;
ub2 = ub;
lb2(j) = ub1(j); 

% printing children nodes
printnode(node.node1, lb1, ub1, depth + 1);
printnode(node.node2, lb2, ub2, depth + 1);
end
