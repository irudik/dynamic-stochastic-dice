%   bounds_calculation: Calculates upper and lower CO2 bounds based off
%   atmospheric CO2 bounds.

% Carbon transfer matrix
c_trans = [81.0712 18.9288 0; 9.7213  85.2787 5; 0 0.3119 99.6881]/100;

% Equations that keep ocean CO2 stocks from jumping up outside the
% collocation domain
c_hi_matrix = [c_trans(2,3) (c_trans(3,3)-1);...
    (c_trans(2,2)-1) c_trans(3,2)];

c_hi_b = [0; -c_trans(1,2)*catm_hi];

% Equations that % Equations that keep ocean CO2 stocks from jumping down
% outside the collocation domain Chi and Clo from jumping down
c_lo_matrix = [c_trans(2,3) (c_trans(3,3)-1);...
     (c_trans(2,2)-1) c_trans(3,2)];

c_lo_b = [0; -c_trans(1,2)*catm_lo];

% Upper bounds for upper and lower oceans
hi23 = c_hi_matrix\c_hi_b;

% Lower bounds for upper and lower oceans
lo23 = c_hi_matrix\c_lo_b;