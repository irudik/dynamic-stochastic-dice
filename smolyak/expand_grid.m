function grid = expand_grid(space)
%   expand_grid: Expands the [-1,1]^d grid into any hyperrectangle
%
%   grid = expand_grid(space);
% INPUTS
%   space: Structure that defines the Smolyak approximation space
% OUTPUTS
%   grid : Grid expanded into the desired hyperrectangle
%
% USES: none
% USED BY: none


center = space.lb+(space.ub-space.lb)/2;
dist = (space.ub-space.lb)/2;
grid = bsxfun(@plus,bsxfun(@times,space.simplex_grid,dist),center);

if size(grid) ~= size(space.simplex_grid)
    error('Grid conformability error.');
end