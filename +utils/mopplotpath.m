%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = mopplotpath(x0, fx0, x1, fx1, psdim, pfdim, opts, varargin)
if opts.PlotPSAxesInd ~= 0 && length(psdim) >= 2
  subplot(opts.PlotAxesPartition1, opts.PlotAxesPartition2, opts.PlotPSAxesInd);
  if length(psdim) == 2
    if isempty(x0)
      plot(x1(psdim(1)), x1(psdim(2)), varargin{:});
    else
      plot([x0(psdim(1)), x1(psdim(1))], [x0(psdim(2)), x1(psdim(2))], varargin{:});
    end
  else
    if isempty(x0)
      plot3(x1(psdim(1)), x1(psdim(2)), x1(psdim(3)), varargin{:});
    else
      plot3([x0(psdim(1)), x1(psdim(1))], [x0(psdim(2)), x1(psdim(2))], [x0(psdim(3)), x1(psdim(3))], varargin{:});
    end
  end
end

if opts.PlotPFAxesInd ~= 0 && length(pfdim) >= 2
  subplot(opts.PlotAxesPartition1, opts.PlotAxesPartition2, opts.PlotPFAxesInd);
  if length(pfdim) == 2
    if isempty(fx0)
      plot(fx1(pfdim(1)), fx1(pfdim(2)), varargin{:});
    else
      plot([fx0(pfdim(1)), fx1(pfdim(1))], [fx0(pfdim(2)), fx1(pfdim(2))], varargin{:});
    end
  else
    if isempty(fx0)
      plot3(fx1(pfdim(1)), fx1(pfdim(2)), fx1(pfdim(3)), varargin{:});
    else
      plot3([fx0(pfdim(1)), fx1(pfdim(1))], [fx0(pfdim(2)), fx1(pfdim(2))], [fx0(pfdim(3)), fx1(pfdim(3))], varargin{:});
    end
  end
end
end