function smol_inds = smol_ind(d,mu)
%   smol_ind: Finds all indices 'i' such that d <= |i=1:d| <= d+mu
%
%   smol_inds = smol_ind(d,mu);
% INPUTS
%   d   : scalar of the number of dimensions
%   mu  : scalar of the approximation level
% OUTPUTS
%   smol_inds: vector of indices that satisfy Smolyak's rule
%
% USES: combinator, poss_inds, uniqueperms
% USED BY: build_grid and poly_ind


% Check for anisotropic grid
if numel(mu) == 1
    max_mu = mu;
elseif length(mu) ~= d
    error('For anisotropic grid, mu must have the same number of elements as the dimension of the grid.');
else
    max_mu = max(mu);
end
    

% Generate initial list of combinations, missing some rows where need to
% flip columns around, will complete in for loop below
poss_inds = combinator(max_mu+1,d,'c','r');

poss_inds = poss_inds((sum(poss_inds,2) <= d+max_mu),:);


poss_inds_tmp = poss_inds;

for i = 1:size(poss_inds,1)
    tmp = uniqueperms(poss_inds_tmp(i,:));
    poss_inds = [poss_inds;tmp];
end

if numel(mu) > 1
    for i = 1:d
        poss_inds(poss_inds(:,i)>mu(i)+1,:) = [];
    end
end
smol_inds = unique(poss_inds, 'rows');