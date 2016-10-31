function [qnew vnew] = maxbell(si,qi,...
    vars,Params,LogicalOpts,space,space_terminal,feedback,weights)

%   maxbell: Initializes state depending variables and maximizes the right
%   hand side of the Bellman.

% Calculate state-dependent exogenous variables for node i
vars.capital = si(1,1)*vars.Ati^(1/(1-Params.kappa))*vars.Lti;
vars.output = vars.Ati*(vars.capital)^Params.kappa*vars.Lti^(1-Params.kappa);
vars.net_output = vars.Ati*(vars.capital)^Params.kappa*vars.Lti^(1-Params.kappa)./(1+Params.b2*si(1,5).^Params.b3);

% Acons*x <= bcons
% Consumption plus abatement cost in output terms less than total output
Acons = [vars.Ati^(1/(1-Params.kappa))*vars.Lti vars.net_output];
bcons = [vars.net_output];

% Lower and upper bounds
% Consumption < output
% Abatement cost < Abatement cost when rate == 100%
xu = [vars.net_output/(vars.Ati^(1/(1-Params.kappa))*vars.Lti); vars.Psii];

% Abatement cost and consumption are weakly positive
xl = [0; 0];


% Set up anonymous functions for maximization, analytic gradient if desired
gradient = @(qi)valuemin_grad(qi,si,vars,Params,LogicalOpts,space,space_terminal,feedback,weights);
value_handle = @(qi)valuemin(qi,si,vars,Params,LogicalOpts,space,space_terminal,feedback,weights,gradient);

% Solver options for if we are in the terminal period
if space.terminal_switch
    options_on = optimoptions(@fmincon,'display', 'off', 'MaxFunEvals', 3000, ...
        'MaxIter', 3000, 'TolFun', 1e-10, 'TolX', 1e-10,'Gradobj', 'off');
    %[qnew, vnew, exitflag, output, lambda, grad, hessian] = fmincon(value_handle,qi',Acons,bcons,[],[],xl,xu,[],options_on);
    
else
    options_on = optimoptions(@fmincon,'display','off', 'MaxFunEvals', 3000, ...
        'MaxIter', 3000, 'TolFun', 1e-6, 'TolX', 1e-10, 'SpecifyObjectiveGradient', false);
    %[qnew, vnew, exitflag, output] = knitromatlab(value_handle,qi,Acons,bcons,[],[],xl,xu,[],[],[],options_file);
    
    %%% To verify gradient is correct un comment this set of options, it
    %%% will check the analytic gradient against a central finite
    %%% difference approximation:
%     options_on = optimoptions(@fmincon,'display','off', 'MaxFunEvals', LogicalOpts.maxfunevals, 'MaxIter', ...
%     LogicalOpts.maxiter, 'TolFun', 1e-6, 'TolX', 1e-10,'SpecifyObjectiveGradient',true,'CheckGradients',true,...
%     'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-7);
   
end

% Minimize the negative right hand side of the Bellman
[qnew, vnew, exitflag, output, lambda, grad, hessian] = ...
    fmincon(value_handle,qi',Acons,bcons,[],[],xl,xu,[],options_on);


qnew=qnew';



end
