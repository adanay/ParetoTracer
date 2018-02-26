% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [multfun] = valmultfun(multfun, doval)
% Validates the multiply functions structure.
% The multiply function must be either a cell array of function handles or 
% a struct. I.e.,
% - multfun = {Jv, wJ, wJv, vH, Hw, Hwv, 
%              Jcv, wJc, wJcv, vHc, Hcw, Hcwv, 
%              Jceqv, wJceq, wJceqv, vHceq, Hceqw, Hceqwv}.
% Note that J stands for Jacobian and H for Hessian.
% Also c stands for inequality constraints and ceq for equality
% constraints.
% - multfun is a struct such that
%   + Jvx = multfun.Jv(x, v) -> J(x) * v
%   + wJx = multfun.wJ(x, w) -> w' * J(x) 
%   + wJvx = multfun.wJv(x, w, v) -> w' * J(x) * v
%   + vHx = multfun.vH(x, v) -> [v' * H1; v' * H2; ...; v' * Hnobj]
%   + Hwx = multfun.Hw(x, w) -> H1 * w1 + H2 * w2 + ... + Hnobj * wnobj
%   + Hwvx = multfun.Hwv(x, w, v) -> (H1 * w1 + H2 * w2 + ... + Hnobj * wnobj) * v
%
%   + Jcvx = multfun.Jcv(x, v) -> Jc(x) * v
%   + wJcx = multfun.wJc(x, w) -> w' * Jc(x)
%   + wJcvx = multfun.wJcv(x, w, v) -> w' * Jc(x) * v
%   + vHcx = multfun.vHc(x, v) -> [v' * Hc1; v' * Hc2; ...; v' * Hcnobj]
%   + Hcwx = multfun.Hcw(x, w) -> Hc1 * w1 + Hc2 * w2 + ... + Hcnobj * wnobj
%   + Hcwvx = multfun.Hcwv(x, w, v) -> (Hc1 * w1 + Hc2 * w2 + ... + Hcnobj * wnobj) * v
%
%   + Jceqvx = multfun.Jceqv(x, v)
%   + wJceqx = multfun.wJceq(x, w)
%   + wJceqvx = multfun.wJceqv(x, w, v)
%   + vHceqx = multfun.vHceq(x, v)
%   + Hceqwx = multfun.Hceqw(x, w)
%   + Hceqwvx = multfun.Hceqwv(x, w, v)

if ~doval 
  % no validation required
  if ~isstruct(multfun)
    multfun = struct('Jv', multfun{1}, 'wJ', multfun{2}, 'wJv', multfun{3}, 'vH', multfun{4}, 'Hw', multfun{5}, 'Hwv', multfun{6},...
                     'Jcv', multfun{7}, 'wJc', multfun{8}, 'wJcv', multfun{9}, 'vHc', multfun{10}, 'Hcw', multfun{11}, 'Hcwv', multfun{12},...
                     'Jceqv', multfun{13}, 'wJceq', multfun{14}, 'wJceqv', multfun{15}, 'vHceq', multfun{16}, 'Hceqw', multfun{17}, 'Hceqwv', multfun{18});
  end
  return
end

if iscell(multfun)
  l = length(multfun);
  
  if l > 0
    Jv = multfun{1};
    if ~isempty(Jv) && ~ischar(Jv) && ~isa(Jv, 'function_handle')
      error('val:valmultfun:InvalidJv', 'The Jacobian multiply function Jv in multfun{1} must be a string or a handle or empty.');
    end
  else
    Jv = [];
  end
  if l > 1
    wJ = multfun{2};
    if ~isempty(wJ) && ~ischar(wJ) && ~isa(wJ, 'function_handle')
      error('val:valmultfun:InvalidwJ', 'The Jacobian multiply function wJ in multfun{2} must be a string or a handle or empty.');
    end
  else
    wJ = [];
  end
  if l > 2
    wJv = multfun{3};
    if ~isempty(wJv) && ~ischar(wJv) && ~isa(wJv, 'function_handle')
      error('val:valmultfun:InvalidwJv', 'The Jacobian multiply function wJv in multfun{3} must be a string or a handle or empty.');
    end
  else
    wJv = [];
  end
  if l > 3
    vH = multfun{4};
    if ~isempty(vH) && ~ischar(vH) && ~isa(vH, 'function_handle')
      error('val:valmultfun:InvalidvH', 'The Hessian multiply function vH in multfun{4} must be a string or a handle or empty.');
    end
  else
    vH = [];
  end
  if l > 4
    Hw = multfun{5};
    if ~isempty(Hw) && ~ischar(Hw) && ~isa(Hw, 'function_handle')
      error('val:valmultfun:InvalidHw', 'The Hessian multiply function Hw in multfun{5} must be a string or a handle or empty.');
    end
  else
    Hw = [];
  end
  if l > 5
    Hwv = multfun{6};
    if ~isempty(Hwv) && ~ischar(Hwv) && ~isa(Hwv, 'function_handle')
      error('val:valmultfun:InvalidHwv', 'The Hessian multiply function Hwv in multfun{6} must be a string or a handle or empty.');
    end
  else
    Hwv = [];
  end
  
  % inequality constraints
  if l > 6
    Jcv = multfun{7};
    if ~isempty(Jcv) && ~ischar(Jcv) && ~isa(Jcv, 'function_handle')
      error('val:valmultfun:InvalidJcv', 'The nonlinear inequality constraints Jacobian multiply function Jcv in multfun{7} must be a string or a handle or empty.');
    end
  else
    Jcv = [];
  end
  if l > 7
    wJc = multfun{8};
    if ~isempty(wJc) && ~ischar(wJc) && ~isa(wJc, 'function_handle')
      error('val:valmultfun:InvalidwJc', 'The nonlinear inequality constraints Jacobian multiply function wJc in multfun{8} must be a string or a handle or empty.');
    end
  else
    wJc = [];
  end
  if l > 8
    wJcv = multfun{9};
    if ~isempty(wJcv) && ~ischar(wJcv) && ~isa(wJcv, 'function_handle')
      error('val:valmultfun:InvalidwJcv', 'The nonlinear inequality constraints Jacobian multiply function wJcv in multfun{9} must be a string or a handle or empty.');
    end
  else
    wJcv = [];
  end
  if l > 9
    vHc = multfun{10};
    if ~isempty(vHc) && ~ischar(vHc) && ~isa(vHc, 'function_handle')
      error('val:valmultfun:InvalidvHc', 'The nonlinear inequality constraints Hessian multiply function vHc in multfun{10} must be a string or a handle or empty.');
    end
  else
    vHc = [];
  end
  if l > 10
    Hcw = multfun{11};
    if ~isempty(Hcw) && ~ischar(Hcw) && ~isa(Hcw, 'function_handle')
      error('val:valmultfun:InvalidHcw', 'The nonlinear inequality constraints Hessian multiply function Hcw in multfun{11} must be a string or a handle or empty.');
    end
  else
    Hcw = [];
  end
  if l > 11
    Hcwv = multfun{12};
    if ~isempty(Hcwv) && ~ischar(Hcwv) && ~isa(Hcwv, 'function_handle')
      error('val:valmultfun:InvalidHcwv', 'The nonlinear inequality constraints Hessian multiply function Hcwv in multfun{12} must be a string or a handle or empty.');
    end
  else
    Hcwv = [];
  end
  
  % equality constraints
  if l > 12
    Jceqv = multfun{13};
    if ~isempty(Jceqv) && ~ischar(Jceqv) && ~isa(Jceqv, 'function_handle')
      error('val:valmultfun:InvalidJceqv', 'The nonlinear equality constraints Jacobian multiply function Jceqv in multfun{13} must be a string or a handle or empty.');
    end
  else
    Jceqv = [];
  end
  if l > 13
    wJceq = multfun{14};
    if ~isempty(wJceq) && ~ischar(wJceq) && ~isa(wJceq, 'function_handle')
      error('val:valmultfun:InvalidwJceq', 'The nonlinear equality constraints Jacobian multiply function wJceq in multfun{14} must be a string or a handle or empty.');
    end
  else
    wJceq = [];
  end
   if l > 14
    wJceqv = multfun{15};
    if ~isempty(wJceqv) && ~ischar(wJceqv) && ~isa(wJceqv, 'function_handle')
      error('val:valmultfun:InvalidwJceqv', 'The nonlinear equality constraints Jacobian multiply function wJceqv in multfun{15} must be a string or a handle or empty.');
    end
  else
    wJceqv = [];
  end
  if l > 15
    vHceq = multfun{16};
    if ~isempty(vHceq) && ~ischar(vHceq) && ~isa(vHceq, 'function_handle')
      error('val:valmultfun:InvalidvHceq', 'The nonlinear equality constraints Hessian multiply function vHceq in multfun{16} must be a string or a handle or empty.');
    end
  else
    vHceq = [];
  end
  if l > 16
    Hceqw = multfun{17};
    if ~isempty(Hceqw) && ~ischar(Hceqw) && ~isa(Hceqw, 'function_handle')
      error('val:valmultfun:InvalidHceqw', 'The nonlinear equality constraints Hessian multiply function Hceqw in multfun{17} must be a string or a handle or empty.');
    end
  else
    Hceqw = [];
  end
  if l > 17
    Hceqwv = multfun{18};
    if ~isempty(Hceqwv) && ~ischar(Hceqwv) && ~isa(Hceqwv, 'function_handle')
      error('val:valmultfun:InvalidHceqwv', 'The nonlinear equality constraints Hessian multiply function Hceqwv in multfun{18} must be a string or a handle or empty.');
    end
  else
    Hceqwv = [];
  end

elseif isstruct(multfun)
  if isfield(multfun, 'Jv')
    Jv = multfun.Jv;
    if ~isempty(Jv) && ~ischar(Jv) && ~isa(Jv, 'function_handle')
      error('val:valmultfun:InvalidJv', 'The Jacobian multiply function multfun.Jv must be a string or a handle or empty.');
    end
  else
    Jv = [];
  end
  if isfield(multfun, 'wJ')
    wJ = multfun.wJ;
    if ~isempty(wJ) && ~ischar(wJ) && ~isa(wJ, 'function_handle')
      error('val:valmultfun:InvalidwJ', 'The Jacobian multiply function multfun.wJ must be a string or a handle or empty.');
    end
  else
    wJ = [];
  end
  if isfield(multfun, 'wJv')
    wJv = multfun.wJv;
    if ~isempty(wJv) && ~ischar(wJv) && ~isa(wJv, 'function_handle')
      error('val:valmultfun:InvalidwJv', 'The Jacobian multiply function multfun.wJv must be a string or a handle or empty.');
    end
  else
    wJv = [];
  end
  if isfield(multfun, 'vH')
    vH = multfun.vH;
    if ~isempty(vH) && ~ischar(vH) && ~isa(vH, 'function_handle')
      error('val:valmultfun:InvalidvH', 'The Hessian multiply function multfun.vH must be a string or a handle or empty.');
    end
  else
    vH = [];
  end
  if isfield(multfun, 'Hw')
    Hw = multfun.Hw;
    if ~isempty(Hw) && ~ischar(Hw) && ~isa(Hw, 'function_handle')
      error('val:valmultfun:InvalidHw', 'The Hessian multiply function multfun.Hw must be a string or a handle or empty.');
    end
  else
    Hw = [];
  end
  if isfield(multfun, 'Hwv')
    Hwv = multfun.Hwv;
    if ~isempty(Hwv) && ~ischar(Hwv) && ~isa(Hwv, 'function_handle')
      error('val:valmultfun:InvalidHwv', 'The Hessian multiply function multfun.Hwv must be a string or a handle or empty.');
    end
  else
    Hwv = [];
  end
  
  % inequality constraints
  if isfield(multfun, 'Jcv')
    Jcv = multfun.Jcv;
    if ~isempty(Jcv) && ~ischar(Jcv) && ~isa(Jcv, 'function_handle')
      error('val:valmultfun:InvalidJcv', 'The nonlinear inequality constraints Jacobian multiply function multfun.Jcv must be a string or a handle or empty.');
    end
  else
    Jcv = [];
  end
  if isfield(multfun, 'wJc')
    wJc = multfun.wJc;
    if ~isempty(wJc) && ~ischar(wJc) && ~isa(wJc, 'function_handle')
      error('val:valmultfun:InvalidwJc', 'The nonlinear inequality constraints Jacobian multiply function multfun.wJc must be a string or a handle or empty.');
    end
  else
    wJc = [];
  end
  if isfield(multfun, 'wJcv')
    wJcv = multfun.wJcv;
    if ~isempty(wJcv) && ~ischar(wJcv) && ~isa(wJcv, 'function_handle')
      error('val:valmultfun:InvalidwJcv', 'The nonlinear inequality constraints Jacobian multiply function multfun.wJcv must be a string or a handle or empty.');
    end
  else
    wJcv = [];
  end
  if isfield(multfun, 'vHc')
    vHc = multfun.vHc;
    if ~isempty(vHc) && ~ischar(vHc) && ~isa(vHc, 'function_handle')
      error('val:valmultfun:InvalidvHc', 'The nonlinear inequality constraints Hessian multiply function multfun.vHc must be a string or a handle or empty.');
    end
  else
    vHc = [];
  end
  if isfield(multfun, 'Hcw')
    Hcw = multfun.Hcw;
    if ~isempty(Hcw) && ~ischar(Hcw) && ~isa(Hcw, 'function_handle')
      error('val:valmultfun:InvalidHcw', 'The nonlinear inequality constraints Hessian multiply function multfun.Hwc must be a string or a handle or empty.');
    end
  else
    Hcw = [];
  end
  if isfield(multfun, 'Hcwv')
    Hcwv = multfun.Hcwv;
    if ~isempty(Hcwv) && ~ischar(Hcwv) && ~isa(Hcwv, 'function_handle')
      error('val:valmultfun:InvalidHcwv', 'The nonlinear inequality constraints Hessian multiply function multfun.Hwcv must be a string or a handle or empty.');
    end
  else
    Hcwv = [];
  end
  
  % equality constraints
  if isfield(multfun, 'Jceqv')
    Jceqv = multfun.Jceqv;
    if ~isempty(Jceqv) && ~ischar(Jceqv) && ~isa(Jceqv, 'function_handle')
      error('val:valmultfun:InvalidJceqv', 'The nonlinear equality constraints Jacobian multiply function multfun.Jceqv must be a string or a handle or empty.');
    end
  else
    Jceqv = [];
  end
  if isfield(multfun, 'wJceq')
    wJceq = multfun.wJceq;
    if ~isempty(wJceq) && ~ischar(wJceq) && ~isa(wJceq, 'function_handle')
      error('val:valmultfun:InvalidwJceq', 'The nonlinear equality constraints Jacobian multiply function multfun.wJceq must be a string or a handle or empty.');
    end
  else
    wJceq = [];
  end
  if isfield(multfun, 'wJceqv')
    wJceqv = multfun.wJceqv;
    if ~isempty(wJceqv) && ~ischar(wJceqv) && ~isa(wJceqv, 'function_handle')
      error('val:valmultfun:InvalidwJceqv', 'The nonlinear equality constraints Jacobian multiply function multfun.wJceqv must be a string or a handle or empty.');
    end
  else
    wJceqv = [];
  end
  if isfield(multfun, 'vHceq')
    vHceq = multfun.vHceq;
    if ~isempty(vHceq) && ~ischar(vHceq) && ~isa(vHceq, 'function_handle')
      error('val:valmultfun:InvalidvHceq', 'The nonlinear equality constraints Hessian multiply function multfun.vHceq must be a string or a handle or empty.');
    end
  else
    vHceq = [];
  end
  if isfield(multfun, 'Hceqw')
    Hceqw = multfun.Hceqw;
    if ~isempty(Hceqw) && ~ischar(Hceqw) && ~isa(Hceqw, 'function_handle')
      error('val:valmultfun:InvalidHceqw', 'The nonlinear equality constraints Hessian multiply function multfun.Hceqw must be a string or a handle or empty.');
    end
  else
    Hceqw = [];
  end
  if isfield(multfun, 'Hceqwv')
    Hceqwv = multfun.Hceqwv;
    if ~isempty(Hceqwv) && ~ischar(Hceqwv) && ~isa(Hceqwv, 'function_handle')
      error('val:valmultfun:InvalidHceqwv', 'The nonlinear equality constraints Hessian multiply function multfun.Hceqwv must be a string or a handle or empty.');
    end
  else
    Hceqwv = [];
  end
else
  if ~isempty(multfun)
    error('val:valmultfun:InvalidMultfun', 'The function multfun must be either a cell array {Jv, wJ, wJv, vH, Hw, Hwv, Jcv, wJc, wJcv, vHc, Hcw, Hcwv, Jceqv, wJceq, wJceqv, vHceq, Hceqw, Hceqwv} or a struct with those fields.');
  end
  Jv = [];
  wJ = [];
  wJv = [];
  vH = [];
  Hw = [];
  Hwv = [];
  Jcv = [];
  wJc = [];
  wJcv = [];
  vHc = [];
  Hcw = [];
  Hcwv = [];
  Jceqv = [];
  wJceq = [];
  wJceqv = [];
  vHceq = [];
  Hceqw = [];
  Hceqwv = [];
end
multfun = struct('Jv', Jv, 'wJ', wJ, 'wJv', wJv, 'vH', vH, 'Hw', Hw, 'Hwv', Hwv,...
                 'Jcv', Jcv, 'wJc', wJc, 'wJcv', wJcv, 'vHc', vHc, 'Hcw', Hcw, 'Hcwv', Hcwv,...
                 'Jceqv', Jceqv, 'wJceq', wJceq, 'wJceqv', wJceqv, 'vHceq', vHceq, 'Hceqw', Hceqw, 'Hceqwv', Hceqwv);
end