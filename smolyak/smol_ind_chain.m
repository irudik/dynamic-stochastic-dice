function phi_chain = smol_ind_chain(n)
%   smol_ind_chain: Generates the disjoint Smolyak indices for the basis
%   functions
%
%   phi_chain = smol_ind_chain(n);
% INPUTS
%   space: Structure that defines the Smolyak approximation space
%   v    : Matrix of values through which to interpolate over the grid
% OUTPUTS
%   coeff: Cell array of disjoint Smolyak indices
%
% USES: none
% USED BY: poly_ind

% Initialize first two
phi_chain{1} = 1;
phi_chain{2} = [2, 3];

curr_val = 4;

for i = 3:n
    end_val = m_i(i);
    temp = curr_val:end_val;
    phi_chain{i} = temp;
    curr_val = end_val+1;
end