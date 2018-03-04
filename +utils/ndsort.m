function [rank, maxrank] = ndsort(pf, s)
% Non-dominated sorting by efficient non-dominated sort (ENS).
%
% rank = ndsort(pf, s) does non-dominated sorting on pf, where pf is the
% matrix of objective values of a set of individuals, and s is the number
% of individuals being sorted at least.
% pf is assumed to be of size (m x nobj) where m is the no of individuals
% and nobj is the no of objectives.
% rank is of size (m x 1) where rank(i) means the front no of the i-th 
% individual. The individuals that have not been sorted are assigned a 
% front no of Inf.
%
% In particular, s = 1 indicates finding only the first non-dominated
% front, s = m / 2 indicates sorting only half of the population (which is 
% often used in the algorithms), and s = Inf indicates sorting the whole 
% population.
%
% [rank, maxrank] = ndsort(pf, s) also returns the maximum front rank 
% besides inf.

[rank, maxrank] = utils.zzz_NDSort(pf, s);
end

