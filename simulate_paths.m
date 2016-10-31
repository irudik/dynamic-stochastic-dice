function [s, qi, v, mus, scc, foc_error, mean_belief, var_belief, channels, fb_out] =...
    simulate_paths(Params,VarsExog,LogicalOpts,space,space_terminal)

% Initialize state vector
s = zeros(Params.timenodes+1,8);
s(1,1) = Params.kstart;
s(1,2) = Params.Mlstart;
s(1,3) = Params.Mustart;
s(1,4) = Params.Mstart;
s(1,5) = Params.Tstart;
s(1,6) = Params.Tostart;
s(1,7) = Params.meanstart;

if LogicalOpts.fb_unc
    s(1,8) = Params.varstart;
else
    s(1,8) = 1e-9;
end

% Start in 2015 if comparing vs Nordhaus
if LogicalOpts.test_nordhaus
    Params.start_time = 2;
    s(2,1) = 172.627/VarsExog.At(2)^(1/(1-Params.kappa))/VarsExog.Lt(2);
    s(2,2) = 18370.477;
    s(2,3) = 1280.635;
    s(2,4) = 863.108;
    s(2,5) = 0.92;
    s(2,6) = 0.043;
    s(2,7) = 0.6;
    s(2,8) = 0;
    Params.fb_cert = 0.6;
    LogicalOpts.fb_unc = 0;
else
    Params.start_time = 1;
end

qstart = [0 0];


fb_out = Params.fb_cert;


% loop through years
for t = Params.start_time:Params.timenodes
    % Set new year
    space.year = (t-1)*Params.timestep;
    
    % Recover time t exogenous variable values
    value_variables = struct(...
            'betaeffi',VarsExog.betaeff,'Bi',VarsExog.Bs(t),...
            'Psii',VarsExog.Psi(t),'eminti',VarsExog.emint(t),...
            'EFi',VarsExog.EF(t),'EFi_next',VarsExog.EF_next(t),...
            'Ati',VarsExog.At(t),'Lti',VarsExog.Lt(t),...
            'Atni',VarsExog.Atn(t),'Ltni',VarsExog.Ltn(t),...
            'sigmai',VarsExog.sigma(t),...
            'eff_adjusti',VarsExog.eff_adjust(t),'eff_adjustni',VarsExog.eff_adjustn(t));
    
    % Get the feedback quadrature at time t
    if LogicalOpts.fb_unc
        
        if s(t,8) > 0
                [feedback_nodes_sim(:,t),feedback_weights_sim(:,t)] = qnwnorm(Params.nqnodes,s(t,7),s(t,8));
        else
            feedback_nodes_sim(:,t) = s(t,7)*ones(Params.nqnodes,1);
            feedback_weights_sim(:,t) = ones(Params.nqnodes,1)/Params.nqnodes;
        end
        
    else
        
        feedback_nodes_sim(:,t) = Params.fb_cert*ones(Params.nqnodes,1);
        feedback_weights_sim(:,t) = ones(Params.nqnodes,1)/Params.nqnodes;
        
    end
    
    % Get 1 dimensional weight vector for our 2 dimensional expectation
    weights_temp = feedback_weights_sim(:,t)*Params.tw';
    weights = reshape(weights_temp,numel(weights_temp),1);

    % Maximize the Bellman
    [qnew, vnew] = maxbell(s(t,:),qstart,...
        value_variables,Params,LogicalOpts,space,space_terminal,feedback_nodes_sim(:,t),weights);
    
    qi(t,:) = qnew;
    v(t,1) = -vnew;
    
    % use this period's control as guess for next period
    qstart = qnew;
    
    % calculate current values of variables that depend on just-chosen
    % optimal actions and revealed parameters
    mu = (qi(t,2)/VarsExog.Psi(t,1))^(1/Params.a2);
    mus(t,1) = mu;
    [snext] = transitions(s(t,:),NaN,qi(t,1),qi(t,2),...
        value_variables,Params,LogicalOpts,space,feedback_nodes_sim(:,t));
    s(t+1,:) = snext;
    forcing = Params.forcing_per_2xco2*(log(s(t+1,4)/Params.Mpre)/log(2))+VarsExog.EF_next(t);
    expected_temp(t,:) = s(t,5) + Params.delay_c1.*(forcing ...
        - (Params.forcing_per_2xco2*((1-feedback_nodes_sim(:,t))./(Params.lambda_0)))*s(t,5) + Params.delay_c3*(s(t,6)-s(t,5)));
    mean_belief(t,:) = feedback_weights_sim(:,t)'*feedback_nodes_sim(:,t);
    var_belief(t,:) = feedback_weights_sim(:,t)'*feedback_nodes_sim(:,t).^2-(feedback_weights_sim(:,t)'*feedback_nodes_sim(:,t)).^2;
end


% converts GtC to ppm
co2conc = s(1:Params.timenodes,2)./Params.gtc_per_ppm;

At = VarsExog.At(1:Params.timenodes,1);
Lt = VarsExog.Lt(1:Params.timenodes,1);

K = s(1:Params.timenodes,1).*At.^(1/(1-Params.kappa)).*Lt;
K_L = s(1:Params.timenodes,1).*At.^(1/(1-Params.kappa));
C = qi(1:Params.timenodes,1).*At.^(1/(1-Params.kappa)).*Lt;
capital = s(1:end-1,1).*At.^(1/(1-Params.kappa)).*Lt;
net_output = At.*(capital).^Params.kappa.*Lt.^(1-Params.kappa)./(1+Params.b2*s(1:end-1,5).^Params.b3);
invest = net_output.*(1-qi(:,2)) - At.^(1/(1-Params.kappa)).*Lt.*qi(:,1);
inv_rate = invest./net_output;

% Put abatement rate in percentage terms
mus = mus*100;

% Calculate SCC and uncertainty channels
if LogicalOpts.scc_switch
    run scc_calculation
else
    scc = 0;
    channels = 0;
    foc_error = 0;
end


% Check if we jump outside the collocation domain during the simulation and
% report what happens
for state = 1:6
    if any(s(1,state) > Params.smax(state))
        warning(['State ' num2str(state) ' jumped upward out of the state space to ' num2str(s(s(:,state)==max(s(:,state)),state))]);
    elseif any(s(:,state) < Params.smin(state))
        warning(['State ' num2str(state) ' jumped downward out of the state space to' num2str(s(s(:,state)==min(s(:,state)),state))]);
    end
end