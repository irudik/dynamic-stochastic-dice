function fct_eval = smol_eval(space,space_terminal,pts,deriv,dist)
%   smol_eval: Evaluates the Smolyak interpolant at a given set of points
%
%   fct_eval = smol_eval(space,pts,deriv);
% INPUTS
%   space: Structure that defines the Smolyak approximation space
%   space_terminal: Structure that defines the Smolyak approximation space
%   for the terminal value function
%   pts  : Points to evaluate the Smolyak interpolant at
%   deriv: Optional switch for evaluating at derivative (doesn't work yet)
% OUTPUTS
%   fct_eval: function evaluated at a certain point
%
% USES: contract_grid
% USED BY: none

if nargin < 4
    if space.terminal_switch
        deriv = zeros(1,space_terminal.d);
    else
        deriv = zeros(1,space.d);
    end
    dist = 0;
end

if size(pts,1) > 1 && size(pts,2) == 1 % for value function
    pts = pts';
end


% If in terminal year, determine if using fixed value or value function
if space.terminal_switch
    if isnumeric(space.terminal_value)
        fct_eval = space.terminal_value*ones(space.nqnodes,1);
    else
        contract_pts = contract_grid(space_terminal,pts);
        contract_B = Smolyak_Polynomial(contract_pts,space_terminal.d,space_terminal.mu_max,space_terminal.pinds,deriv,dist);
        
        % Function value
        fct_eval = contract_B*space_terminal.c_store;
        
        % Warn if put in something illegal
        if ~strcmp(space.terminal_value,'VAL')
            warning('Terminal year value is neither a fixed value nor a value function, defaulting to a value function.');
        end
    end
else
    
    % Contract given points to unit simplex
    contract_pts = contract_grid(space,pts);
    
    % Create B matrix of polynomials at given points
    contract_B = Smolyak_Polynomial(contract_pts,space.d,space.mu_max,space.pinds,deriv,dist);

    % Function value
    fct_eval = contract_B*space.c_store{space.year/space.timestep+2};
end