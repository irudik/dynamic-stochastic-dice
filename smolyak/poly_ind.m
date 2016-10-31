function poly_inds = poly_ind(d, mu, inds)
%   poly_ind: Computes the indices for all tensor products of Chebyshev
%   polynomials to interpolate using Smolyak method
%
%   poly_inds = poly_ind(d, mu, inds);
% INPUTS
%   d        : Scalar of the number of dimensions
%   mu       : Scalar of the approximation level
%   inds     : Optional vector of indices that satisfy Smolyak's rule
% OUTPUTS
%   poly_inds: Indices of tensor products to use for interpolation
%
% USES: smol_inds, smol_ind_chain, smol_chain, cartprodt
% USED BY: cheb_matrix

% If passed empty, calculate here
if isempty(inds)
    inds = smol_ind(d,mu);
end

% Check for anisotropic
if numel(mu) == 1
    max_mu = mu;
else
    max_mu = max(mu);
end

smol_chain = smol_ind_chain(max_mu+1);
poly_inds = [];

% Loop down satisfactory smolov indices
for i = 1:size(inds,1)
    clear temp
    
    %Loop across index sum, pull values out of chain for tensor product
    for j = 1:size(inds,2)
        temp{j} = smol_chain{inds(i,j)};
    end
    
    poly_inds = [poly_inds;cartprodt(temp)];
end

poly_inds = unique(poly_inds,'rows');