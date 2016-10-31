function  [c_out,q_out] = value_approx(topdir,approximation_level,...
    final_year,terminal_value,terminal_path,LogicalOpts,Params)

%   value_approx: Builds the collocation grid. Initializes the model
%   parameters and exogenous variable trajectories. Steps backwards through
%   time from the terminal period. Runs simulations.

close all;

startcount = tic;

% Write information to a .txt file on where collocation points jumped
% outside the domain
LogicalOpts.outbounds = 1;        % =1 write to .txt
LogicalOpts.outboundswarn = 0;  % print warning to screen every time jump outside state space  =1 yes, = 0 no
outboundsfile = 'domain_out';    % prefix of file to which to write info

% Build parameter structure
run('parameterization1.m');

%%  Begin initialization

% set working directory
if ~exist(topdir, 'dir')
    mkdir(topdir);
end
cd(topdir);
disp(['topdir: ',topdir]);

% Make the directory for the results folder
dirname_tail = [' unc ' num2str(LogicalOpts.fb_unc) ' learn ' num2str(LogicalOpts.fb_learn) ...
     ' approx_level ' num2str(approximation_level)];
dirname = [Params.test_prefix dirname_tail];
disp(['directory: \',dirname]);

% Directory for preliminary model information
cd(topdir);
preldir=['prelim ' dirname];
if ~exist(preldir, 'dir')
    mkdir(preldir);
end
cd(preldir);
save('prelim_workspace');


disp(['Beginning iteration process.'])

% Build approximation space
tic
run build_approx_space
s = space.grid;
toc

% Pre-compute exogenous variable trajectories
VarsExog = transitions_exogenous(Params,LogicalOpts);

% Build feedback quadrature at each node
% Must build a matrix of quadrature nodes to accomodate 0 variance states
for node = 1:size(s,1)
    
    if LogicalOpts.fb_unc
        
        % If non-zero variance use quadrature packages
        if s(node,8) > 0
            
            [feedback_nodes(:,node),feedback_weights(:,node)] = qnwnorm(Params.nqnodes,s(node,7),s(node,8));
            
        else
            
            % If 0 variance, just plug in mean
            feedback_nodes(:,node) = s(node,7)*ones(Params.nqnodes,1);
            feedback_weights(:,node) = ones(Params.nqnodes,1)/Params.nqnodes;
            
        end
        
    else
        
        % If perfect info use the mean from the truncated quadrature
        feedback_nodes(:,node) = s(node,7)*ones(Params.nqnodes,1);
        feedback_weights(:,node) = ones(Params.nqnodes,1)/Params.nqnodes;
        
    end
    
    
end

% Step backwards through the finite horizon problem
for period=Params.final_year/Params.timestep:-1:0
    
    % Set the correct year for the exogenous trajectories and for
    % coefficients of the value function in space.c_store
    space.year = period*Params.timestep;
    
    % initialization of vector of optimal values
    v = zeros(nodes_total,1);
    
    % Switch for if we are using the terminal value function, turn off once
    % we are before the terminal period
    if space.year ~= Params.final_year
        space.terminal_switch = 0;
    end
    
    % If we are in the terminal year
    if space.terminal_switch
        
        if ~isnumeric(space.terminal_value)
            % Load terminal value function
            load(terminal_path);
        else
            % Fixed payoff
            space_terminal = 0;
        end
        
    end
    
    % Loop over collocation grid points and maximize the Bellman
    [q,v,c] = ...
        function_iterate([topdir,'\',preldir,'\',outboundsfile,'_'],...
        s,v,q,VarsExog,LogicalOpts,Params,space,space_terminal,feedback_nodes,feedback_weights,period);
    
    % Store coefficients
    space.c_store{period+1,1} = c;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STORAGE AND PLOTTING %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(topdir)

dirname = [dirname];

if ~exist(dirname, 'dir')
    mkdir(dirname)
end
cd(dirname)

% move outbounds files to final directory
try
    movefile([topdir,'\',preldir,'\',outboundsfile,'*']);
    clear wrotefiles;
catch
end


time_elapsed=toc(startcount);
disp(['Time to complete: ',num2str(time_elapsed/60),' minutes.'])

save('workspace')

tic;
% Will write total time to file
time_total = sum(time_elapsed(1:end))/60;

% return final coefficients and approximation space
c_out = c;
q_out = q;


%% Plotting

if LogicalOpts.plot==1
    
    cd(topdir)
    if exist([dirname '\workspace.mat'], 'file')==2
        if exist(preldir, 'dir')
            [status,message,messid] = rmdir(preldir,'s');
            disp([messid ' ' preldir])
            clear status message messid
        end
    end
    
    disp(['Setting up to plot.']);
    cd(topdir)
    cd(dirname)
    run time_paths
    cd(topdir)
    cd(dirname)
    
else % don't plot
    
    toc;
    cd(topdir)
    if exist([dirname '\workspace.mat'], 'file')==2
        if exist(preldir, 'dir')
            [status,message,messid] = rmdir(preldir,'s');
            disp([messid ' ' preldir])
            clear status message messid
        end
    end
    time_total = sum(time_elapsed(1:end))/60;
    disp(['It took ',num2str(time_total),' minutes to reach the stored basis function coefficients,']);
    
end

cd(topdir)
cd(dirname)

end