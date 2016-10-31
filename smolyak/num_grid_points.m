function num_points = num_grid_points(d,mu)
%   num_grid_points: Returns the number of grid points for a given
%   dimension and approximation level up to mu=3
%
%   num_points = num_grid_points(d,mu);
% INPUTS
%   d         : scalar of the number of dimensions
%   mu        : scalar of the approximation level
% OUTPUTS
%   num_points: Number of grid points
%
% USES: none
% USED BY: none

switch mu
    case 1
        num_points = 2*d+1;
    case 2
        num_points = 1+4*d+4*d*(d-1)/2;
    case 3
        num_points = 1 + 8*d + 12*d*(d-1)/2. + 8*d*(d-1)*(d-2)/6;
end