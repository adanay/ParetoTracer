% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [multfun] = valhessmultfun(multfun, doval)
% Validates the Hessian multiply functions structure.
% The multiply function must be either a cell array of function handles or 
% a struct. I.e.,
% - multfun = {vH, Hw, Hwv, 
%              vHc, Hcw, Hcwv, 
%              vHceq, Hceqw, Hceqwv}.
% - multfun is a struct such that
%   + vHx = multfun.vH(x, v) => [v' * H1; v' * H2; ...; v' * Hnobj]
%   + Hwx = multfun.Hw(x, w) => H1 * w1 + H2 * w2 + ... + Hnobj * wnobj
%   + Hwvx = multfun.Hwv(x, w, v) => (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v
%   + vHcx = multfun.vHc(x, v) => [v' * Hc1; v' * Hc2; ...; v' * Hcnobj]
%   + Hcwx = multfun.Hcw(x, w) => Hc1 * w1 + Hc2 * w2 + ... + Hcnobj * wnobj
%   + Hcwvx = multfun.Hcwv(x, w, v) => (Hc1 * w1 + Hc2 * w2 + ... + Hcnobj * wnobj) * v
%   + vHceqx = multfun.vHceq(x, v)
%   + Hceqwx = multfun.Hceqw(x, w)
%   + Hceqwvx = multfun.Hceqwv(x, w, v)

if ~doval 
  % no validation required
  if ~isstruct(multfun)
    multfun = struct('vH', multfun{1}, 'Hw', multfun{2}, 'Hwv', multfun{3},...
                     'vHc', multfun{4}, 'Hcw', multfun{5}, 'Hcwv', multfun{6},...
                     'vHceq', multfun{7}, 'Hceqw', multfun{8}, 'Hceqwv', multfun{9});
  end
  return
end

if iscell(multfun)
  l = length(multfun);
  
  if l > 0
    vH = multfun{1};
    if ~isempty(vH) && ~ischar(vH) && ~isa(vH, 'function_handle')
      error('val:valhessmultfun:InvalidvH', 'The Hessian multiply function vH in multfun{1} must be a string or a handle or empty.');
    end
  else
    vH = [];
  end
  if l > 1
    Hw = multfun{2};
    if ~isempty(Hw) && ~ischar(Hw) && ~isa(Hw, 'function_handle')
      error('val:valhessmultfun:InvalidHw', 'The Hessian multiply function Hw in multfun{2} must be a string or a handle or empty.');
    end
  else
    Hw = [];
  end
  if l > 2
    Hwv = multfun{3};
    if ~isempty(Hwv) && ~ischar(Hwv) && ~isa(Hwv, 'function_handle')
      error('val:valhessmultfun:InvalidHwv', 'The Hessian multiply function Hwv in multfun{3} must be a string or a handle or empty.');
    end
  else
    Hwv = [];
  end
  
  % inequality constraints
  if l > 3
    vHc = multfun{4};
    if ~isempty(vHc) && ~ischar(vHc) && ~isa(vHc, 'function_handle')
      error('val:valhessmultfun:InvalidvHc', 'The nonlinear inequality constraints Hessian multiply function vHc in multfun{4} must be a string or a handle or empty.');
    end
  else
    vHc = [];
  end
  if l > 4
    Hcw = multfun{5};
    if ~isempty(Hcw) && ~ischar(Hcw) && ~isa(Hcw, 'function_handle')
      error('val:valhessmultfun:InvalidHcw', 'The nonlinear inequality constraints Hessian multiply function Hcw in multfun{5} must be a string or a handle or empty.');
    end
  else
    Hcw = [];
  end
  if l > 5
    Hcwv = multfun{6};
    if ~isempty(Hcwv) && ~ischar(Hcwv) && ~isa(Hcwv, 'function_handle')
      error('val:valhessmultfun:InvalidHcwv', 'The nonlinear inequality constraints Hessian multiply function Hcwv in multfun{6} must be a string or a handle or empty.');
    end
  else
    Hcwv = [];
  end
  
  % equality constraints
  if l > 6
    vHceq = multfun{7};
    if ~isempty(vHceq) && ~ischar(vHceq) && ~isa(vHceq, 'function_handle')
      error('val:valhessmultfun:InvalidvHceq', 'The nonlinear equality constraints Hessian multiply function vHceq in multfun{7} must be a string or a handle or empty.');
    end
  else
    vHceq = [];
  end
  if l > 7
    Hceqw = multfun{8};
    if ~isempty(Hceqw) && ~ischar(Hceqw) && ~isa(Hceqw, 'function_handle')
      error('val:valhessmultfun:InvalidHceqw', 'The nonlinear equality constraints Hessian multiply function Hceqw in multfun{8} must be a string or a handle or empty.');
    end
  else
    Hceqw = [];
  end
  if l > 8
    Hceqwv = multfun{9};
    if ~isempty(Hceqwv) && ~ischar(Hceqwv) && ~isa(Hceqwv, 'function_handle')
      error('val:valhessmultfun:InvalidHceqwv', 'The nonlinear equality constraints Hessian multiply function Hceqwv in multfun{9} must be a string or a handle or empty.');
    end
  else
    Hceqwv = [];
  end

elseif isstruct(multfun)
  if isfield(multfun, 'vH')
    vH = multfun.vH;
    if ~isempty(vH) && ~ischar(vH) && ~isa(vH, 'function_handle')
      error('val:valhessmultfun:InvalidvH', 'The Hessian multiply function multfun.vH must be a string or a handle or empty.');
    end
  else
    vH = [];
  end
  if isfield(multfun, 'Hw')
    Hw = multfun.Hw;
    if ~isempty(Hw) && ~ischar(Hw) && ~isa(Hw, 'function_handle')
      error('val:valhessmultfun:InvalidHw', 'The Hessian multiply function multfun.Hw must be a string or a handle or empty.');
    end
  else
    Hw = [];
  end
  if isfield(multfun, 'Hwv')
    Hwv = multfun.Hwv;
    if ~isempty(Hwv) && ~ischar(Hwv) && ~isa(Hwv, 'function_handle')
      error('val:valhessmultfun:InvalidHwv', 'The Hessian multiply function multfun.Hwv must be a string or a handle or empty.');
    end
  else
    Hwv = [];
  end
  
  % inequality constraints
  if isfield(multfun, 'vHc')
    vHc = multfun.vHc;
    if ~isempty(vHc) && ~ischar(vHc) && ~isa(vHc, 'function_handle')
      error('val:valhessmultfun:InvalidvHc', 'The nonlinear inequality constraints Hessian multiply function multfun.vHc must be a string or a handle or empty.');
    end
  else
    vHc = [];
  end
  if isfield(multfun, 'Hcw')
    Hcw = multfun.Hcw;
    if ~isempty(Hcw) && ~ischar(Hcw) && ~isa(Hcw, 'function_handle')
      error('val:valhessmultfun:InvalidHcw', 'The nonlinear inequality constraints Hessian multiply function multfun.Hwc must be a string or a handle or empty.');
    end
  else
    Hcw = [];
  end
  if isfield(multfun, 'Hcwv')
    Hcwv = multfun.Hcwv;
    if ~isempty(Hcwv) && ~ischar(Hcwv) && ~isa(Hcwv, 'function_handle')
      error('val:valhessmultfun:InvalidHcwv', 'The nonlinear inequality constraints Hessian multiply function multfun.Hwcv must be a string or a handle or empty.');
    end
  else
    Hcwv = [];
  end
  
  % equality constraints
  if isfield(multfun, 'vHceq')
    vHceq = multfun.vHceq;
    if ~isempty(vHceq) && ~ischar(vHceq) && ~isa(vHceq, 'function_handle')
      error('val:valhessmultfun:InvalidvHceq', 'The nonlinear equality constraints Hessian multiply function multfun.vHceq must be a string or a handle or empty.');
    end
  else
    vHceq = [];
  end
  if isfield(multfun, 'Hceqw')
    Hceqw = multfun.Hceqw;
    if ~isempty(Hceqw) && ~ischar(Hceqw) && ~isa(Hceqw, 'function_handle')
      error('val:valhessmultfun:InvalidHceqw', 'The nonlinear equality constraints Hessian multiply function multfun.Hceqw must be a string or a handle or empty.');
    end
  else
    Hceqw = [];
  end
  if isfield(multfun, 'Hceqwv')
    Hceqwv = multfun.Hceqwv;
    if ~isempty(Hceqwv) && ~ischar(Hceqwv) && ~isa(Hceqwv, 'function_handle')
      error('val:valhessmultfun:InvalidHceqwv', 'The nonlinear equality constraints Hessian multiply function multfun.Hceqwv must be a string or a handle or empty.');
    end
  else
    Hceqwv = [];
  end
else
  if ~isempty(multfun)
    error('val:valhessmultfun:InvalidMultfun', 'The function multfun must be either a cell array {vH, Hw, Hwv, Hcw, vHc, Hcwv, vHceq, Hceqw, Hceqwv} or a struct with those fields.');
  end
  vH = [];
  Hw = [];
  Hwv = [];
  vHc = [];
  Hcw = [];
  Hcwv = [];
  vHceq = [];
  Hceqw = [];
  Hceqwv = [];
end
multfun = struct('vH', vH, 'Hw', Hw, 'Hwv', Hwv,...
                 'vHc', vHc, 'Hcw', Hcw, 'Hcwv', Hcwv,...
                 'vHceq', vHceq, 'Hceqw', Hceqw, 'Hceqwv', Hceqwv);
end