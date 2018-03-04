function [] = plotpareto(ps, pf, opts, varargin)
% Plots the Pareto set and front.
% ps is the Pareto set.
% pf is the Pareto front.
% The opts structure have the following fields:
% - PlotAxesPartition1, PlotAxesPartition2: Breaks the Figure window into 
%   a m-by-n matrix of small axes. By default it is (1 x 2).
% - PlotPSAxesInd, PlotPFAxesInd: Indices of the axes where the PS and the  
%   PF will be plotted. 0 means that they are not plotted. By default they
%   are 2 and 1 respectively.
% - PlotPSDims, PlotPFDims: Respectively the dimensions of the PS and PF to
%   be plotted. By default they are both 1 : 3.
% It can be empty.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

n = size(ps, 2);
nobj = size(pf, 2);

if ~exist('opts', 'var') || isempty(opts)
  opts = struct(...
    'PlotAxesPartition1', 1,...
    'PlotAxesPartition2', 2,...
    'PlotPSAxesInd', 2,...% 0 means do not plot
    'PlotPFAxesInd', 1,...% 0 means do not plot
    'PlotPSDims', 1 : min(n, 3),...
    'PlotPFDims', 1 : min(nobj, 3)); 
end

if isfield(opts, 'PlotPSDims')
  psdim = opts.PlotPSDims;
  psdim(psdim > n | psdim < 1) = [];
  psdim = psdim(1 : min([end, 3]));
else
  psdim = 1 : min(n, 3);
end
ps = ps(:, psdim);
n = size(ps, 2);

if isfield(opts, 'PlotPFDims')
  pfdim = opts.PlotPFDims;
  pfdim(pfdim > nobj | pfdim < 1) = [];
  pfdim = pfdim(1 : min([end, 3]));
else
  pfdim = 1 : min(nobj, 3);
end
pf = pf(:, pfdim);
nobj = size(pf, 2);

if isfield(opts, 'PSAxesInd')
  psaxesind = opts.PSAxesInd;
else
  psaxesind = 2;
end

if isfield(opts, 'PFAxesInd')
  pfaxesind = opts.PlotPFAxesInd;
else
  pfaxesind = 1;
end

if isfield(opts, 'PlotAxesPartition1')
  axespart1 = opts.PlotAxesPartition1;
else
  axespart1 = 1;
end

if isfield(opts, 'PlotAxesPartition2')
  axespart2 = opts.PlotAxesPartition2;
else
  if psaxesind ~= 0 && pfaxesind ~= 0
    axespart2 = 2;
  else
    axespart2 = 1;
  end
end

% variable side
if psaxesind ~= 0
  if n < 2
    error('utils:plotpareto:NLessThan2', 'There must be at least 2 variables to plot.');
  end
  
  plotparetoset();
end

% objective side
if pfaxesind ~= 0
  if nobj < 2
    error('utils:plotpareto:NObjLessThan2', 'There must be at least 2 objectives to plot.');
  end
  
  plotparetofront();
end

hold off

  function [] = plotparetoset()
    subplot(axespart1, axespart2, psaxesind);
    hold on
    axis square
    xlabel(['$x_{', num2str(psdim(1)), '}$'], 'fontsize', 18, 'interpreter','latex');
    ylabel(['$x_{', num2str(psdim(2)), '}$'], 'fontsize', 18, 'interpreter','latex');
    if length(psdim) == 2
      plot(ps(:, 1), ps(:, 2), varargin{:});
    else
      zlabel(['$x_{', num2str(psdim(3)), '}$'], 'fontsize', 18, 'interpreter','latex');
      view(3);
      grid on
      plot3(ps(:, 1), ps(:, 2), ps(:, 3), varargin{:});
    end
  end

  function [] = plotparetofront()
    subplot(axespart1, axespart2, pfaxesind);
    hold on
    axis square
    xlabel(['$f_{', num2str(pfdim(1)), '}$'], 'fontsize', 18, 'interpreter','latex');
    ylabel(['$f_{', num2str(pfdim(2)), '}$'], 'fontsize', 18, 'interpreter','latex');
    if length(pfdim) == 2
      plot(pf(:, 1), pf(:, 2), varargin{:});
    else
      zlabel(['$f_{', num2str(pfdim(3)), '}$'], 'fontsize', 18, 'interpreter','latex');
      view(3);
      grid on
      plot3(pf(:, 1), pf(:, 2), pf(:, 3), varargin{:});
    end
  end
end



