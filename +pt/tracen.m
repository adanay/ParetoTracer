function [itlast, result, stats, EXITFLAG] = tracen(it0, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, opts, stats, pcoutfcn)
% Performs the continuation for MOPs with 3 or more objectives.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if ~exist('pcoutfcn', 'var')
  pcoutfcn = [];
end

diropts = 'IncNegDir';

i = 0; % queue first index
j = 1; % queue last index

% queue initial setup
N = opts.InitSetSize / 4;
pending0(N) = pt.traceit();
pending1(N) = pt.traceit();

pending0(1) = it0;
pending1(1) = it1;

while true
  i = i + 1; % index of the starting point 
  
  it0 = pending0(i);
  it1 = pending1(i);
  
  % predictors -> correctors
  [it1, itcs, result, stats, stop] = pt.predictcorrect(it0, it1, result, objfun, x0, funvals0, lb, ub, lincon, nonlcon, multfun, diropts, opts, stats, pcoutfcn);
  
  % number of correctors
  m = length(itcs);
  
  % last point
  if m == 0
    itlast = it1;
  else
    itlast = itcs(end);
  end
  
  % stop condition (output function)
  if stop
    EXITFLAG = -1;
    return
  end
  
  % reset
  for k = 1 : m
    % resizing sets
    if j == N
      N = ceil(N * opts.SetSizeGrowFact + 1);
      pending0(1 : i - 1) = [];
      pending1(1 : i - 1) = [];
      pending0(N) = pt.traceit();
      pending1(N) = pt.traceit();
      j = j - i + 1;
      i = 1;
    end
    
    j = j + 1;
    
    pending0(j) = it1;
    pending1(j) = itcs(k);
  end
  
  % stop condition (no more points found)
  if i >= j
    EXITFLAG = 1;
    return
  end
  
  % stop condition (max number of iterations reached)
  if stats.PCIts == opts.PCMaxIts
    EXITFLAG = 0;
    return
  end
end
end
