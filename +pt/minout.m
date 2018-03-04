%

% Copyright (c) 2018 Adanay Mart�n & Oliver Sch�tze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = minout(info)

if info.opts.SuppressOutput || info.opts.OptSuppressOutput
  return
end

if ~info.opts.OptSuppressPrint
  pt.minprint(info)
end

if ~info.opts.OptSuppressPlot  
  pt.minplot(info)
end
end

