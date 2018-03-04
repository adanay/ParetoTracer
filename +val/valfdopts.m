function [validopts] = valfdopts(opts, n, doval)
% Validates the FD options.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if ~doval
  % no validation required
  validopts = opts;
  return
end

validopts = fd.defopts();
if isempty(opts)
  return
end

if ~isstruct(opts)
  error('val:valfdopts:NonStructOpts', 'The options (opts) must be a structure.'); 
end

vfn = fieldnames(validopts);
fn = fieldnames(opts);

for i = 1 : length(vfn)
  index = find(strcmpi(fn, vfn{i})); 
  if ~isempty(index) % user value was specified
    value = opts.(fn{index(1)}); % user value
    if ~isempty(value)  
      if val.valfdoptnv(vfn{i}, value, n) % validate user value
        validopts.(vfn{i}) = value; % assign user value
      else
        error('val:valfdopts:InvalidOptValue', ' The option %s has an invalid value.', vfn{i});
      end
    end
  end
end

if validopts.FDMinChange > validopts.FDMaxChange 
  error('val:valfdopts:FDMinMaxChange', 'The option FDMinChange is greater than the option FDMaxChange.');
end
end

