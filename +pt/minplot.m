%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = minplot(info)

if strcmpi(info.opts.OptPlotMode, 'off')
  return
end

switch (upper(info.PHASEFLAG))
  case 'INIT'
    minplotinit(info);
  case 'FINIT'
    minplotfinit(info);
  otherwise
    minplotit(info);
end
end

function [] = minplotinit(info)
[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(info.it1);
sizes = {n, nobj, na, naeq, nc, nceq};
N = n + nobj + na + naeq + nc + nceq + n + n;

psdim = info.opts.PlotPSDims;
psdim(psdim > N | psdim < 1) = [];
psdim = psdim(1 : min([end, 3]));

maxdim = max(psdim);

pfdim = info.opts.PlotPFDims;
pfdim(pfdim > nobj | pfdim < 1) = [];
pfdim = pfdim(1 : min([end, 3]));

% variables
if info.opts.PlotPSAxesInd ~= 0 && length(psdim) >= 2
  subplot(info.opts.PlotAxesPartition1, info.opts.PlotAxesPartition2, info.opts.PlotPSAxesInd);
  hold on
  axis square
  if length(psdim) == 2
    box on
  else
    grid on
    view(3);
  end
  utils.mopvarlabels(psdim, sizes);
end

% objectives
if info.opts.PlotPFAxesInd ~= 0 && length(pfdim) >= 2
  subplot(info.opts.PlotAxesPartition1, info.opts.PlotAxesPartition2, info.opts.PlotPFAxesInd);
  hold on
  axis square
  if length(pfdim) == 2
    box on
  else
    view(3);
    grid on
  end
  utils.mopobjlabels(pfdim);
end

% first point
if ~strcmpi(info.opts.OptPlotMode, 'iter')
  return
end

args = {'marker', '^', 'color', [0, 0.4470, 0.7410]};
argsit = {'marker', '.', 'color', 'black', 'markersize', 8};

x1 = utils.mopaugmentx(info.it1.x, info.it1.w, maxdim);
fx1 = info.it1.fx;

utils.mopplotpath([], [], x1, fx1, psdim, pfdim, info.opts, args{:});
utils.mopplotpath([], [], x1, fx1, psdim, pfdim, info.opts, argsit{:});
drawnow
end

function [] = minplotfinit(info)
[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(info.result);
N = n + nobj + na + naeq + nc + nceq + n + n;

psdim = info.opts.PlotPSDims;
psdim(psdim > N | psdim < 1) = [];
psdim = psdim(1 : min([end, 3]));

maxdim = max(psdim);

pfdim = info.opts.PlotPFDims;
pfdim(pfdim > nobj | pfdim < 1) = [];
pfdim = pfdim(1 : min([end, 3]));

% last point
args = {'marker', 'o', 'color', [0.6350, 0.0780, 0.1840]};
argsit = {'marker', '.', 'color', 'black', 'markersize', 8};

x = utils.mopaugmentx(info.result.x, info.result.w, maxdim);
fx = info.result.fx;

utils.mopplotpath([], [], x, fx, psdim, pfdim, info.opts, args{:});
utils.mopplotpath([], [], x, fx, psdim, pfdim, info.opts, argsit{:});
drawnow
end

function [] = minplotit(info)
if ~strcmpi(info.opts.OptPlotMode, 'iter')
  return
end
if ~strcmpi(info.PHASEFLAG, 'STEP')
  return
end

[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(info.it1);
N = n + nobj + na + naeq + nc + nceq + n + n;

psdim = info.opts.PlotPSDims;
psdim(psdim > N | psdim < 1) = [];
psdim = psdim(1 : min([end, 3]));

maxdim = max(psdim);

pfdim = info.opts.PlotPFDims;
pfdim(pfdim > nobj | pfdim < 1) = [];
pfdim = pfdim(1 : min([end, 3]));

args = {'marker', '.', 'color', 'black', 'markersize', 8};

x0 = utils.mopaugmentx(info.it1.x, info.it1.w, maxdim);
fx0 = info.it1.fx;

x1 = utils.mopaugmentx(info.it2.x, info.it2.w, maxdim);
fx1 = info.it2.fx;

utils.mopplotpath(x0, fx0, x1, fx1, psdim, pfdim, info.opts, args{:});
drawnow
end

