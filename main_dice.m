%   main_dice: The main script for the model. It contains the main
%   switches and parameter values that are varied in the paper.
%
% This model requires use of the CompEcon toolbox by Miranda and Fackler.
% It can be downloaded at: http://www4.ncsu.edu/~pfackler/compecon/toolbox.html

% Clear workspace
clear all;
clear global;
close all;

% Define root directory
dir_root = fileparts(mfilename('fullpath'));
addpath(genpath([dir_root]));
addpath(genpath([dir_root '/smolyak']));
addpath(genpath([dir_root '/compecon']));
addpath(genpath([dir_root '/results']));

% Leading name of directory
Params.test_prefix = 'TESTING';

%%% Feedback parameter options %%%

% Switch for feedback uncertainty
LogicalOpts.fb_unc = 1;

% Switch for feedback learning
LogicalOpts.fb_learn = 0;

% Upper and lower bounds on the mean of the climate sensitivty distribution
Params.cs_upper = 6;
Params.cs_lower = 2;

%%% Approximation options %%%

% Smolyak approximation level
approximation_level = 2*ones(1,8);

% Number of Gauss-Hermite quadrature nodes (total number is Params.nqnodes^2)
Params.nqnodes = 7;

% Preference parameters
Params.rho = 1 - (2); % 1 - 1/EIS
Params.risk = 1 - (2); % 1 - RRA

% Final year of horizon, must end in 5 since the initial year is 2005 and
% the timestep is 10 years
final_year = 2555;

% Scrap value function after final_year
terminal_value = 'VAL'; % enter a number if constant continuation value at terminal time
                        % enter 'VAL' if using value function from infinite horizon problem
                        
% Path to terminal value function
terminal_path = [dir_root '\terminal_workspaces\terminal_2_550.mat'];

% Move to root directory
cd(dir_root);

% Plot outcomes
LogicalOpts.plot = 1;

%%%%%%% END PARAMETER DEFINITIONS  %%%%%%%%%%
initialize_model(dir_root,approximation_level,final_year,terminal_value,terminal_path,LogicalOpts,Params);