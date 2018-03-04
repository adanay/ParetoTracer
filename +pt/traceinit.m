function [it0, it1, result, opts, stats] = traceinit(objfun, lincon, nonlcon, x0, funvals0, opts)
% Initializes the structures required by the Pareto Tracer algorithm. 

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% statistics
stats = pt.tracestats();
    
% previous iteration values
it0 = pt.traceit();

% current iteration values
it1 = it0;
it1.x = x0;

% available values
if ~isempty(funvals0)
  it1.fx = funvals0.fx;
  it1.Jx = funvals0.Jx;
  it1.Hx = funvals0.Hx;
  it1.Hident = funvals0.Hident;
  it1.Lx = funvals0.Lx;
  it1.Lmodif = funvals0.Lmodif;
  it1.ax = funvals0.ax;
  it1.aeqx = funvals0.aeqx;
  it1.cx = funvals0.cx;
  it1.ceqx = funvals0.ceqx;
  it1.Jcx = funvals0.Jcx;
  it1.Hcx = funvals0.Hcx;
  it1.Jceqx = funvals0.Jceqx;
  it1.Hceqx = funvals0.Hceqx;
  it1.dceqx = funvals0.dceqx;
  it1.dcx = funvals0.dcx;
  it1.w = funvals0.w;
  it1.wfx = funvals0.fx;
  it1.wJx = funvals0.Jx;
  it1.Hwx = funvals0.Hx;
  it1.wax = funvals0.ax;
  it1.waeqx = funvals0.aeqx;
  it1.wcx = funvals0.cx;
  it1.wceqx = funvals0.ceqx;
  it1.wJcx = funvals0.Jcx;
  it1.Hcwx = funvals0.Hcx;
  it1.wJceqx = funvals0.Jceqx;
  it1.Hceqwx = funvals0.Hceqx;
end

% ensures function evaluations are computed
[it1, stats] = pt.feval(it1, objfun, lincon, nonlcon, false, opts, stats);

% initializes the result structure
[it1, result, opts, stats] = pt.traceresult(it1, opts, stats);
end

