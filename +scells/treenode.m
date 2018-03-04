classdef treenode < handle
  % Defines a node class for a tree.
  
  % Copyright (c) 2018 Adanay Martín & Oliver Schütze.
  % This file is subject to the terms and conditions defined in
  % the file 'LICENSE.txt', which is part of this source code package.

  properties
    node1 % left node
    node2 % right node
    index % index in an external list (optional)
  end
end
