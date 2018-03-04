function [mu, R] = mud(Mu, J, d)
% Computes the mu vectors given a set of vectors in objective space.
% A linear system of equations is solved for this purpose:
%
% [-J * Mu] * mu = [d] (ensures that the mu components sum up to 0)    
% [1 ... 1]        [0]
%
% d is of size (m x nobj) where nobj is the number of objectives and m is
% the number of directions we want to follow in objective space.

% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

% sizes
nobj = size(J, 1);
N = size(Mu, 2);
m = size(d, 1);

W = [J * Mu; [ones(1, nobj), zeros(1, N - nobj)]];
D = [d'; zeros(1, m)];

% W = J * Mu;
% D = d';

o.RECT = true;
[mu, R] = linsolve(W, D, o);
% if R < nobj % this is not correct: R may be the cond no
%   warning('pt:mud:rankDeficientMatrix', 'The matrix utilized to compute mu (i.e., W = [J * Mu; [1...1 0...0]) is rank deficient with rank=%d and nobj=%d.', R, nobj);
% end
mu = mu';
end


