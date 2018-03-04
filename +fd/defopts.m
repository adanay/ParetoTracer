function [opts] = defopts()
% Default options for finite differences (FD). 
% - FDType: FD, used to estimate gradients, are either: 
%   + 'forward' (the default), or 
%   + 'central', which takes twice as many function evaluations but should 
%     be more accurate.
%   Note: Central FD is not currently supported.
% - TypicalX: Scalar or vector that specifies typical magnitude of  
%   variables. The default is 1.
% - FDMinChange: Minimum change allowed in variables. Default to 0.
% - FDMaxChange: Maximum change allowed in variables. Default to Inf.
% - FDStepSize: Scalar or vector step size factor. When you  set it to a 
%   vector t, the change in variables is calculated as: 
%   e = t .* max(abs(x), abs(TypicalX)) .* sign'(x); 
%   where sign'(x) is -1 for x < 0 and 1 for x >= 0. 
%   The default is sqrt(eps).
% - FDStepSize2: Same as above but used for central finite differences or 
%   when the Jacobian is not provided when approximated Hessians. For 
%   central finite differences, when set it to a vector t, the change in 
%   variables is calculated as:
%   e = t .* max(abs(x), abs(TypicalX));
%   The default is eps^(1/3). 
% - UseVectorized: Defines whether multiple values of the function can be
%   obtained with a single call to f, e.g., Y = f(X) where X is a matrix 
%   where each row represents a single individual, and analogously, Y is a
%   matrix where each row represents a function evaluation. For the  
%   Jacobian function it is similar. If X has a size of (m x n), Y = J(X)  
%   will have a size of (m x nobj x n). The default is false.
% - LargeScale: If the problem is large scale, a matrix of size (n x n)
%   will never be formed, just the resulting (nobj x n) Jacobian. 
%   When approximating Hessians this setting applies only to directional 
%   derivatives. If the problem is large scale, a matrix of size (n x n) 
%   will never be formed, just the resulting (nobj x n) directional 
%   derivative. The default is false.  
% - FunValCheck: Checks whether the objective function values are valid. If 
%   true, displays an error when the objective function returns a value  
%   that is complex, NaN, or Inf. 
%   When approximating Hessians, if the Jacobian is not valid, a warning 
%   will be thrown and the Hessian will be approximated to the identity.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

opts = struct(...
  'FDType', 'forward',...
  'TypicalX', 1,...
  'FDMinChange', 0,...
  'FDMaxChange', Inf,...
  'FDStepSize', sqrt(eps),...
  'FDStepSize2', eps^(1/3),... % used for central FD or to approximate Hessians from the objective function
  'UseVectorized', false,...
  'LargeScale', false,...
  'FunValCheck', true,...
  'ValidateInput', true);
end
