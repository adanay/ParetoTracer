%

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = mopobjlabels(pfdim)
l = length(pfdim);

xlabel(['$f_{', num2str(pfdim(1)), '}$'], 'fontsize', 18, 'interpreter','latex');
if l > 1
  ylabel(['$f_{', num2str(pfdim(2)), '}$'], 'fontsize', 18, 'interpreter','latex');
end
if l > 2
  zlabel(['$f_{', num2str(pfdim(3)), '}$'], 'fontsize', 18, 'interpreter','latex');
end
end