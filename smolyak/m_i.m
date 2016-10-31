function poly_deg = m_i(i)
%   m_i: Returns number of points for 'i'th nested set
%
%   poly_deg = m_i(i);
% INPUTS
%   i       : Index of the nested set
% OUTPUTS
%   poly_deg: Number of points in the indexed set
%
% USES: none
% USED BY: cheb_matrix, smol_ind_chain


if i < 0
    error('Need non-negative integer i.');
elseif i == 0
    poly_deg = 0;
elseif i == 1
    poly_deg = 1;
else
    poly_deg = 2^(i-1)+1;
end