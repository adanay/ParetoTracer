%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = traceout(info)

if info.opts.SuppressOutput || info.opts.PCSuppressOutput
  return
end

if ~info.opts.PCSuppressPrint
  pt.traceprint(info)
end

if ~info.opts.PCSuppressPlot  
  pt.traceplot(info)
end
end

