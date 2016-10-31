function extrema_set = s_n(n)
%   s_n: Returns the set S_n which is the n'th set of Chebyshev extrema
%
%   extrema_set = s_n(n);
% INPUTS
%   n          : Index of the desired set of Chebyshev extrema
% OUTPUTS
%   extrema_set: Vector of Chebyshev extrema

if n==1
    extrema_set = 0;
end

% if extrema

M_i = 2^(n-1)+1; % may not need the +1 for root nodes
comp_vals = 1:M_i;

% else 
%   M_i = 2^(n-1);
%   comp_vals = 1:M_i;
% end

extrema_set = -1*cos(pi*(comp_vals-1)/(M_i-1));
extrema_set(abs(extrema_set) < 1e-14) = 0;