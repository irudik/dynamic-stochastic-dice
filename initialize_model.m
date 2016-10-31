function initialize_model(dir_root,approximation_level,...
    final_year,terminal_value,terminal_path,LogicalOpts,Params)
%   initialize_model: Initializes the domain of the model and sets up
%   directories.


% save results in
dir_main = dir_root;
dir_root = [dir_root,'\results'];

% Translate climate sensitivity into feedback
Params.lambda_0 = 1.2;
Params.fb_upper = 1-Params.lambda_0/Params.cs_upper;
Params.fb_lower = 1-Params.lambda_0/Params.cs_lower;

% Order of states is:
% Effective capital, lower ocean CO2, upper ocean CO2, atmospheric CO2,
% surface temperature, ocean temperature, mean and variance of beliefs


% Given selections for the atmospheric CO2 domain, use the
% transition equations to infer what the domain for the ocean CO2
% stocks must be so that don't transition outside the domain in those
% states
catm_hi = 1700;
catm_lo = 580;
run bounds_calculation

% Some nodes very slightly leave the domain due to numerical error.
% Adjust the ocean CO2 domains so they are the ones where this
% happens and not atmospheric CO2. This has no impact on the
% results.
lo23 = fliplr(lo23'+[100 100]);
hi23 = fliplr(hi23'-[1000 6000]);


Params.smin = [1.70 lo23 catm_lo 00.00 00.00 Params.fb_lower 0.000];
Params.smax = [6.00 hi23 catm_hi 10.60 10.60 Params.fb_upper 0.13^2];

% Display domain
disp(['State space: min ' mat2str(Params.smin) ', max ' mat2str(Params.smax)]);

% Set up directories
topdir = [dir_root];

% Proceed to value function approximation
[c,q] = value_approx(topdir,approximation_level,final_year,terminal_value,...
    terminal_path,LogicalOpts,Params);
