function grid = contract_grid(space,pts)
%   contract_grid: Projects approximated function evaluation points back
%   into the unit simplex. Each row is a different point. Each column is a
%   different dimension.
%
%   grid = contract_grid(space,pts);
% INPUTS
%   space: Structure that defines the Smolyak approximation space
%   pts  : Bounds of the grid to contract into [-1,1]^d
% OUTPUTS
%   grid : Contracted grid in [-1,1]^d
%
% USES: none
% USED BY: smol_eval

% Contract points. Catch exception if potentially simulating using an
% approximation before adaptive grid was implemented.

grid = bsxfun(@rdivide,bsxfun(@minus,pts,space.center),space.dist);