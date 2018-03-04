function [validopts] = valptopts(opts, n, doval)
% Validates the Pareto Tracer options.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

if ~doval
  % no validation required
  validopts = opts;
  return
end

validopts = pt.defopts();

if isempty(opts)
  return
end

if ~isstruct(opts)
  error('val:valptopts:NonStructOpts', 'The options (opts) must be a structure.'); 
end

% valid options regarding FD
validfdopts = val.valfdopts(opts, n, doval);
vfdfn = fieldnames(validfdopts);

vfn = fieldnames(validopts);
fn = fieldnames(opts);

for i = 1 : length(vfn)
  fdindex = find(strcmpi(vfdfn, vfn{i}));
  if ~isempty(fdindex) % FD option
    value = validfdopts.(vfdfn{fdindex(1)}); % FD value
    validopts.(vfn{i}) = value; % assign FD value
  else
    index = find(strcmpi(fn, vfn{i}));
    if ~isempty(index) % user value was specified
      value = opts.(fn{index(1)}); % user value
      if ~isempty(value)
        % validate user value
        if val.valptoptnv(vfn{i}, value, n)
          validopts.(vfn{i}) = value; % assign user value
        else
          error('val:valptopts:InvalidOptValue', ' The option %s has an invalid value.', vfn{i});
        end
      end
    end
  end
end
end

