%   build_approx_space: This is the upper level script to build the Smolyak
%   grid and interpolant. This script builds the approximation structure
%   'space.'. It contains all things necessary for updating coefficients
%   (besides new function evaluations) and evaluating the interpolant at
%   arbitrary points. Some lower level scripts of this Smolyak
%   approximation method is a MATLAB adaptation of Python code for the Judd
%   et al. (2014) version of the Smolyak method by Chase Coleman and
%   Spencer Lyon at https://github.com/EconForge. The main scripts to
%   construct the Smolyak polynomals is an altered version of the code
%   provided by Lilia and Serguei Maliar. The Smolyak_Polynomial.m file is
%   edited to increase speed and allow for analytic (non-mixed) derivatives.
%   Note: combinator, cartprodt, and uniqueperms were obtained from the
%   MathWorks file exchange.
%
% INPUTS-------------------------------------------------------------------
%   approximation_level: Desired approximation level
%
%   Params.smax        : Vector of upper bounds of the approximation space
%   dimensions
%
%   Params.smax        : Vector of lower bounds of the approximation space
%   dimensions
%
% OUTPUTS------------------------------------------------------------------
%   space              : Structure containing the grid and pre-computed
%   inverted matrix of Chebyshev polynomials evaluated at the Smolyak grid
%   points in order to update the interpolant coefficients

mu = approximation_level;


d = length(Params.smin);
space.lb = Params.smin;
space.ub = Params.smax;


% Error check
if d < 1
    error('Must have non-negative integer dimensions.');
elseif mu < 1
    error('Must have non-negative integer approximation level.');
end


% Store Smolyak grid features in structure
space.d = d;
space.mu = mu;
space.mu_max = max(mu);
space.inds = smol_ind(d,mu);
space.pinds = poly_ind(d,mu,space.inds);
space.num_nodes = size(space.pinds,1);


% Determine center and radius of approx space
space.center = space.lb+(space.ub-space.lb)/2;
space.dist = (space.ub-space.lb)/2;

% Generate simplex grid
space.iso_inds = Smolyak_Elem_Isotrop(space.d,space.mu_max);
space.pinds = Smolyak_Elem_Anisotrop(space.iso_inds,space.mu);
space.simplex_grid = Smolyak_Grid(d,space.mu_max,space.pinds);

% Generate 'B' matrix of Cheb polynomials evaluated at all grid points
space.B = Smolyak_Polynomial(space.simplex_grid,d,space.mu_max,space.pinds);

% Set elements equal to zero that are slightly off due to numerical
% error
space.B(abs(space.B)<1e-10) = 0;

% Expand to approximation space
space.grid = expand_grid(space);

% Store grid in state matrix
disp(['Grid constructed with ' num2str(size(space.grid,1)) ' nodes.']);
coeffs_total = space.num_nodes;
nodes_total = space.num_nodes;

% Pre-compute inverse to save time for future coefficient calculations
space.Binv = space.B\eye(size(space.B));

% Initializes control guesses
try
    q = space_terminal.q;
catch
    q = [0.5*ones(nodes_total,1) .001*ones(nodes_total,1)];
end

% Store key variables in space structure
space.final_year = Params.final_year;
space.terminal_switch = 1;
space.terminal_value = terminal_value;
space.terminal_path = terminal_path;
space.timestep = Params.timestep;