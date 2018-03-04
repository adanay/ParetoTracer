%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = traceplot(info)

if strcmpi(info.opts.PCPlotMode, 'off')
  return
end

switch (upper(info.PHASEFLAG))
  case 'INIT'
    traceplotinit(info);
  case 'FINIT'
    traceplotfinit(info);
  otherwise
    traceplotit(info);
end
end

function traceplotinit(info)
if ~(info.opts.OptSuppressOutput || info.opts.OptSuppressPlot || strcmpi(info.opts.OptPlotMode, 'off'))
  return
end

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
if ~strcmpi(info.opts.PCPlotMode, 'iter') || maxdim > n
  args = {'marker', '*', 'color', 'black', 'MarkerEdgeColor', [0, 0.4470, 0.7410]};
else
  args = {'marker', '^', 'color', 'black', 'MarkerEdgeColor', [0, 0.4470, 0.7410]};
end

x1 = utils.mopaugmentx(info.it1.x, info.it1.w, maxdim);
fx1 = info.it1.fx;
utils.mopplotpath([], [], x1, fx1, psdim, pfdim, info.opts, args{:});

if strcmpi(info.opts.PCPlotMode, 'iter') && maxdim <= n
  args = {'marker', '.', 'color', 'black', 'markersize', 8};
  utils.mopplotpath([], [], x1, fx1, psdim, pfdim, info.opts, args{:});
end

drawnow
end

function traceplotfinit(info)
[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(info.it1);
N = n + nobj + na + naeq + nc + nceq + n + n;

psdim = info.opts.PlotPSDims;
psdim(psdim > N | psdim < 1) = [];
psdim = psdim(1 : min([end, 3]));

maxdim = max(psdim);

pfdim = info.opts.PlotPFDims;
pfdim(pfdim > nobj | pfdim < 1) = [];
pfdim = pfdim(1 : min([end, 3]));

% last point
if ~strcmpi(info.opts.PCPlotMode, 'iter') || maxdim > n
  args = {'marker', '*', 'color', [0.6350, 0.0780, 0.1840]};
else
  args = {'marker', 'o', 'color', [0.6350, 0.0780, 0.1840]};
end
x1 = utils.mopaugmentx(info.it1.x, info.it1.w, maxdim);
fx1 = info.it1.fx;
utils.mopplotpath([], [], x1, fx1, psdim, pfdim, info.opts, args{:});
drawnow
end

function [] = traceplotit(info)
[n, nobj, na, naeq, nc, nceq] = pt.mopsizes(info.it1);
N = n + nobj + na + naeq + nc + nceq + n + n;

psdim = info.opts.PlotPSDims;
psdim(psdim > N | psdim < 1) = [];
psdim = psdim(1 : min([end, 3]));

maxdim = max(psdim);

pfdim = info.opts.PlotPFDims;
pfdim(pfdim > nobj | pfdim < 1) = [];
pfdim = pfdim(1 : min([end, 3]));

switch (upper(info.PHASEFLAG))
  case 'PRED'
    traceplotpred();
  case 'CORRECT'
    traceplotcorrect();
  case 'SAVE'
    if info.opts.PCDrawCells && ~isempty(info.result.tree)
      tracedrawcell();
    end
end

  function [] = traceplotpred()
    if ~strcmpi(info.opts.PCPlotMode, 'iter') || maxdim > n
      return
    end
    
    args = {'color', [0.4940, 0.1840, 0.5560]};
    
    x0 = info.it1.x;
    fx0 = info.it1.fx;
    
    x1 = info.itp.x;
    fx1 = info.itp.fx;
    
    utils.mopplotpath(x0, fx0, x1, fx1, psdim, pfdim, info.opts, args{:});
    drawnow
  end

  function [] = traceplotcorrect()
    if isempty(info.it0.x) && info.stats.PCIts == 0 && ~(info.opts.OptSuppressOutput || info.opts.OptSuppressPlot || strcmpi(info.opts.OptPlotMode, 'off'))
      return
    end 
      
    % first point
    if isempty(info.it0.x) && info.stats.PCIts == 0
      if ~strcmpi(info.opts.PCPlotMode, 'iter') || maxdim > n
        args = {'marker', '*', 'color', 'black', 'MarkerEdgeColor', [0, 0.4470, 0.7410]};
      else
        args = {'marker', '^', 'color', 'black', 'MarkerEdgeColor', [0, 0.4470, 0.7410]};
      end
      
      x0 = utils.mopaugmentx(info.itp.x, info.it1.w, maxdim);
      fx0 = info.itp.fx;

      x1 = utils.mopaugmentx(info.itc.x, info.it1.w, maxdim);
      fx1 = info.itc.fx;
      utils.mopplotpath(x0, fx0, x1, fx1, psdim, pfdim, info.opts, args{:});
    end
    
    if strcmpi(info.opts.PCPlotMode, 'iter') && maxdim <= n
      if isempty(info.it0.x) && info.stats.PCIts == 0
        args = {'color', 'black'};
      else
        args = {'color', [0.4660, 0.6740, 0.1880]}; 
      end
      
      x0 = info.itp.x;
      fx0 = info.itp.fx;
      
      x1 = info.itc.x;
      fx1 = info.itc.fx;
      
      utils.mopplotpath(x0, fx0, x1, fx1, psdim, pfdim, info.opts, args{:});
      
      args = {'marker', '.', 'color', 'black', 'markersize', 8};
    else
      args = {'marker', info.opts.PCMarker, 'color', 'black'};
    end
    
    if strcmpi(info.opts.PCPlotMode, 'flow')
      x0 = utils.mopaugmentx(info.it1.x, info.it1.w, maxdim);
      fx0 = info.it1.fx;
    else
      x0 = [];
      fx0 = [];
    end
    
    x1 = utils.mopaugmentx(info.itc.x, info.itc.w, maxdim);
    fx1 = info.itc.fx;
    
    utils.mopplotpath(x0, fx0, x1, fx1, psdim, pfdim, info.opts, args{:});
    drawnow
  end

  function [] = tracedrawcell()
    if info.opts.PlotPFAxesInd ~= 0
      subplot(info.opts.PlotAxesPartition1, info.opts.PlotAxesPartition2, info.opts.PlotPFAxesInd);
      hold on
      c = info.itc.cell(pfdim);
      r = info.result.cellradius(pfdim);
      if info.itc.cellinserted
        scells.drawcell(c, r, 'facealpha', 0.2, 'facecolor', 'y');
        drawnow
      end
    end
  end
end




