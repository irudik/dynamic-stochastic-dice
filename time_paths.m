% time_paths: simulates the trajectories for a model

% Define user as 'ivan' or 'derek'
user = 'ivan';

% Define root directory
dir_root = fileparts(mfilename('fullpath'));
addpath(genpath([dir_root]));
addpath(genpath([dir_root '/smolyak']));
addpath(genpath([dir_root '/compecon']));
addpath(genpath([dir_root '/results']));
topdir = [dir_root '\results'];

% Test accuracy against Nordhaus paths (begin in 2015)
LogicalOpts.test_nordhaus = 0;

% Calculate SCC
LogicalOpts.scc_switch = 0;

% Terminal simulation period
T = 200;

% Initialize vectors, etc
run initialize_simulation
Params.final_year = T;

% Pre-allocate exogenous variable trajectories
[VarsExog] = ...
    transitions_exogenous(Params,LogicalOpts);

% Parallel simulation of the trajectories
% Simulate trajectories
[s_out, q_out, v_out, mus_out, scc_out, ...
    foc_error_out, mean_belief_out, var_belief_out, channels_out, fb_out]...
    = simulate_paths(Params,VarsExog,LogicalOpts,...
    space,space_terminal);

k_store = s_out(:,1);
mlo_store = s_out(:,2);
mup_store = s_out(:,3);
matm_store = s_out(:,4);
ts_store = s_out(:,5);
to_store = s_out(:,6);
mean_store = s_out(:,7);
var_store = s_out(:,8);
q_store = q_out;
v_store = v_out;
mus_store = mus_out;
scc_store = scc_out;
foc_store = foc_error_out;
mean_belief_store = mean_belief_out;
var_belief_store = var_belief_out;
fb_store = fb_out;


save_string = ['nord_' num2str(LogicalOpts.test_nordhaus) '.mat'];
cd([topdir '\' dirname '\time_paths']);

clear logical_*
save(save_string);

clear options x1 qstart value_handle temp_adder

% change abatement to percentage

time_timepath = toc(time_path_count);
time_timepath = time_timepath/60;
disp(['Time to complete plotting time paths: ',num2str(time_timepath),' minutes'])
