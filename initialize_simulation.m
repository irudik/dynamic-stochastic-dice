%   initialize_simulation: initializes the simulation, makes folder to
%   store simulated trajectories.

clear s q

% Start timer
time_path_count = tic;

% Initial count for jumping out during simulation, may do it because of
% noise or bad approximation domain
jumpcount=0;

% Make simulation directory
if ~exist([topdir '\' dirname '\time_paths'])
    mkdir([topdir '\' dirname '\time_paths']);
end
cd([topdir '\' dirname '\time_paths']);

% Times at which we will simulate the problem
ts = (0:Params.timestep:T)';

% The number of periods we are simulating
Params.timenodes = size(ts,1);

% Initial guess for utility maximization (i.e., in period 1)
qstart=[0 0];