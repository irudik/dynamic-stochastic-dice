function coeff = smol_fit(space,v)
%   smol_fit: Solves for the coefficients of the Smolyak interpolant. This
%   is the value function update step of value function iteration.
%
%   coeff = smol_fit(space,v);
% INPUTS
%   space: Structure that defines the Smolyak approximation space
%   v    : Matrix of values through which to interpolate over the grid
% OUTPUTS
%   coeff: Coefficients of the Smolyak interpolant for values v
%
% USES: none
% USED BY: none


if size(space.Binv,1) ~= size(space.Binv,2) || size(space.Binv,1) ~= length(v)
    error('Conformation error. B matrix isnt square or doesnt match length of evaluated grid.');
end

if size(v,1) == 1
    v = v';
end

coeff = space.Binv*v;