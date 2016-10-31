% scc_calculation: Calculates the optimal carbon tax and decomposes it into
% the different channels described in the paper.

% Generate quadrature of next period's temperature along simulated path %
forcing = Params.forcing_per_2xco2*(log(s(2:end,4)/Params.Mpre)/log(2))+ VarsExog.EF_next;
clear temp_next exp_next

% Derivative vectors for continuation value
k_deriv = [1 0 0 0 0 0 0 0];
m_deriv = [0 0 0 1 0 0 0 0];
t_deriv = [0 0 0 0 1 0 0 0];
t_deriv_2 = [0 0 0 0 2 0 0 0];
t_deriv_3 = [0 0 0 0 3 0 0 0];
mu_deriv = [0 0 0 0 0 0 1 0];
mu_deriv_2 = [0 0 0 0 0 0 2 0];
var_deriv = [0 0 0 0 0 0 0 1];

% Generate cs nodes and next period's expected temperature for SCC, one
% period ahead
% calculation % # nodes x # periods
for i = 1:size(s,1)-1
    
    % Beliefs about next period's temperature. Forcing is already shifted
    % forward one period.
    temp_nodes(:,i) = (s(i,5) + Params.delay_c1.*(forcing(i) ...
        - (Params.forcing_per_2xco2.*((1-feedback_nodes_sim(:,i))./(Params.lambda_0))).*s(i,5) + Params.delay_c3*(s(i,6)-s(i,5))))';
    
    % All combinations of temperature and shock quadrature for double
    % expectation
    temp_comb = combvec(temp_nodes(:,i)',Params.te')';
    
    % Quadrature for next period's temperature (double expectation)
    temp_next(:,i) = temp_comb(:,1) + temp_comb(:,2);

    
    if LogicalOpts.fb_learn
        
        H = temp_next(:,i) - (s(i,5) + Params.delay_c1.*(forcing(i) ...
            - (Params.forcing_per_2xco2./(Params.lambda_0))*s(i,5) + Params.delay_c3*(s(i,6)-s(i,5))));
        
        gamma_param(i,1) = (Params.delay_c1*Params.forcing_per_2xco2*s(i,5)/Params.lambda_0);
        
        % Next period's feedback mean
        exp_next(:,i) = (s(i,8)*(gamma_param(i))*(H) + Params.tvar*s(i,7))/...
            (Params.tvar+s(i,8)*(gamma_param(i)).^2);
       
        % Change in variance from an increase in temperature, for the
        % active learning channel
        dsigma_dt(:,i) = -2*s(i,8)^2*Params.tvar*(gamma_param(i))*(Params.delay_c1*Params.forcing_per_2xco2/Params.lambda_0)/(Params.tvar+s(i,8)*(gamma_param(i)).^2)^2;

        
    else
        exp_next(:,i) = Params.meanstart*ones(size(temp_next(:,i)));
        gamma_param(i,1) = (Params.delay_c1*Params.forcing_per_2xco2*s(i,5)/Params.lambda_0);
    end
end



% Marginal effect of emissions today on temperature tomorrow,
% deterministics in DICE for one period ahead (t+2 is uncertain)
% # nodes x # periods
dT_de = repmat((Params.forcing_per_2xco2/log(2))*Params.delay_c1./s(2:Params.timenodes,4)',size(temp_next,1),1);

% Expected temperature and mean of feedback distribution
exp_td(:,1) = weights'*temp_next;
exp_ed(:,1) = weights'*exp_next; % should be equal to current belief

% Variance of temperature and mean of feedback distribution
var_td(:,1) = weights'*(bsxfun(@minus,temp_next, weights'*temp_next)).^2;
var_ed(:,1) = weights'*((bsxfun(@minus,exp_next, weights'*exp_next)).^2);

% Covariance between temperature and mean of feedback distribution
cov_t_e(:,1) = weights'*(bsxfun(@minus,exp_next, weights'*exp_next).*bsxfun(@minus,temp_next, weights'*temp_next));

% Next period's states for the value function channels
for i = 2:size(s,1)-1
    
    % Time t+1 state quadrature
    s_base{i-1} = [(bsxfun(@times,s(i,1:4),ones(size(temp_next,1),4))) ...
        temp_next(:,i-1) bsxfun(@times,s(i,6),ones(size(temp_next,1),1)) ...
        exp_next(:,i-1) bsxfun(@times,s(i,8),ones(size(temp_next,1),1))];
    
    % Expected state
    s_base_exp{i-1} = [s(i,1:4) exp_td(i-1) s(i,6) exp_ed(i-1) s(i,8)];
    
    % Certainty-equivalent state, time t prior mean, zero
    % variance
    s_base_ce{i-1} = [s(i,1:4) exp_td(i-1) s(i,6) s(i-1,7) 0];
end


% Two period ahead calculations, requires maximizing the Bellman at our
% quadrature nodes for time t+1 states and then getting the quadrature for
% time t+2 states, all are uncertain since time t+1 controls are uncertain
% In principle, all channels could be calculated this way, but we are only
% using it for active learning
for t = 2:size(s_base,2)+1 % # loop through years
    
    % Initialize year and exogenous trajectories
    space.year = (t-1)*Params.timestep;
    index_valueparams = t;
    run value_params
    
    parfor i = 1:size(s_base{1},1) % # Loop over t+1 state quadrature nodes
        
        % Build quadrature for feedback, s_base{t} is the set quadrature
        % nodes of time t+1's state with the time t information set
        feedback_nodes_quad = qnwnorm(Params.nqnodes,s_base{t-1}(i,7),s_base{t-1}(i,8));
        
        % Maximize the Bellman at this quadrature node for the state in
        % time t+1
        [q_out, v_out] = maxbell(s_base{t-1}(i,:),[0 0],...
            value_variables,Params,LogicalOpts,space,space_terminal,feedback_nodes_quad,weights);
        mu_out = (q_out(2)/VarsExog.Psi(t,1))^(1/Params.a2);
        
        % Get 2 period ahead state conditional on the quadrature node for
        % time t+1
        [snext] = transitions(s_base{t-1}(i,:),NaN,q_out(1),q_out(2),...
        value_variables,Params,LogicalOpts,space,feedback_nodes_quad);
    
        % Store states that transition from t+1 to t+2 deterministically
        m_quad(:,i,t) = snext(4)*ones(Params.nqnodes^2,1);
        mup_quad(:,i,t) = snext(3)*ones(Params.nqnodes^2,1);
        mlo_quad(:,i,t) = snext(2)*ones(Params.nqnodes^2,1);
        k_quad(:,i,t) = snext(1)*ones(Params.nqnodes^2,1);
        var_quad(:,i,t) = snext(8)*ones(Params.nqnodes^2,1);
        to_quad(:,i,t) = snext(5)*ones(Params.nqnodes^2,1);

        % Mean belief and temperature require another expectation since
        % transition to t+1 and then to t+2 are both uncertain
        
        % Forcing at time t+2
        forcing_next = Params.forcing_per_2xco2*(log((snext(4))/Params.Mpre)/log(2))+value_variables.EFi_next;
        
        % Quadrature of t+2 temperature given t+1 quadrature node we're at
        temp_next = (s_base{t-1}(i,5) + ...
            Params.delay_c1.*(forcing_next - (Params.forcing_per_2xco2*(1-feedback_nodes_quad)/(Params.lambda_0))*s_base{t-1}(i,5) + ...
            Params.delay_c3*(s_base{t-1}(i,6)-s_base{t-1}(i,5))))';
        
        % dT_t+2/dT_t+1
        dT_dT = 1-Params.delay_c1*Params.forcing_per_2xco2/Params.lambda_0*(1-feedback_nodes_quad) + Params.delay_c1*Params.delay_c3;
        
        temp_comb = combvec(temp_next,Params.te')';
        
        % Quadrature of t+2 temperature given t+1 quadrature node we're at
        temp_quad(:,i,t) = temp_comb(:,1) + temp_comb(:,2);
        
        temp_comb = combvec(dT_dT',Params.te')';
        
        dT_dT_2(:,i,t) = temp_comb(:,1);
        
        % Non-random part of temperature transition
        H = temp_quad(:,i,t) - (s_base{t-1}(i,5) + Params.delay_c1.*(forcing_next ...
            - (Params.forcing_per_2xco2./(Params.lambda_0))*s_base{t-1}(i,5) + Params.delay_c3*(s_base{t-1}(i,6)-s_base{t-1}(i,5))));
        H_store(:,i,t) = H;
        % Constants multiplying the uncertain feedback term
        gamma2 = (Params.delay_c1*Params.forcing_per_2xco2*s_base{t-1}(i,5)/Params.lambda_0);
        
        % Quadrature of time t+2 mean given t+1 quadrature node, should
        % equal the initial belief used in qnwnorm if we took the
        % expectation over this, and then time t belief if we take the
        % expectation over time t+1 quadrature nodes
        exp_quad(:,i,t) = (s_base{t-1}(i,8)*(gamma2)*(H) + Params.tvar*s_base{t-1}(i,7))/...
            (Params.tvar+s_base{t-1}(i,8)*(gamma2).^2);
        
        if ~LogicalOpts.fb_learn
            exp_quad(:,i,t) = s_base{t-1}(i,7);
        end
    end
    
end

% Remove zero row, need to do this for parallelization in MATLAB above
m_quad(:,:,1) = [];
k_quad(:,:,1) = [];
mlo_quad(:,:,1) = [];
mup_quad(:,:,1) = [];
var_quad(:,:,1) = [];
exp_quad(:,:,1) = [];
to_quad(:,:,1) = [];
temp_quad(:,:,1) = [];
dT_dT_2(:,:,1) = [];
H_store(:,:,1) = [];

% Expectations two periods ahead, we can use the same weights since they
% are the same each period, but the location of the quadrature points
% change over time which is captured in the 3 dimensional array
for t = 1:size(exp_quad,3)
    exp_m_2(t,1) = weights'*m_quad(:,:,t)*weights;
    exp_mup_2(t,1) = weights'*mup_quad(:,:,t)*weights;
    exp_mlo_2(t,1) = weights'*mlo_quad(:,:,t)*weights;
    exp_to_2(t,1) = weights'*to_quad(:,:,t)*weights;
    exp_var_2(t,1) = weights'*var_quad(:,:,t)*weights;
    exp_k_2(t,1) = weights'*k_quad(:,:,t)*weights;
    exp_exp_2(t,1) = weights'*exp_quad(:,:,t)*weights;
    exp_temp_2(t,1) = weights'*temp_quad(:,:,t)*weights;
end

% Taking 4 expectations now so get new weighting vector
weights2 = weights*weights';
weights2 = reshape(weights2,numel(weights2),1);

% Put in two dimensional form, reshaping the two quadrature dimensions of
% the state at time t+2 into one quadrature dimension so we can plug it
% into a state matrix for value function evaluation
temp_quad = reshape(temp_quad,size(temp_quad,1)*size(temp_quad,2),size(temp_quad,3))';
k_quad = reshape(k_quad,size(m_quad,1)*size(m_quad,2),size(m_quad,3))';
mup_quad = reshape(mup_quad,size(m_quad,1)*size(m_quad,2),size(m_quad,3))';
mlo_quad = reshape(mlo_quad,size(m_quad,1)*size(m_quad,2),size(m_quad,3))';
to_quad = reshape(to_quad,size(m_quad,1)*size(m_quad,2),size(m_quad,3))';
var_quad = reshape(var_quad,size(m_quad,1)*size(m_quad,2),size(m_quad,3))';
exp_quad = reshape(exp_quad,size(exp_quad,1)*size(exp_quad,2),size(exp_quad,3))';
m_quad = reshape(m_quad,size(m_quad,1)*size(m_quad,2),size(m_quad,3))';
dT_dT_2 = reshape(dT_dT_2,size(dT_dT_2,1)*size(dT_dT_2,2),size(dT_dT_2,3))';
H_quad = reshape(H_store,size(H_store,1)*size(H_store,2),size(H_store,3))';

% For the mean active learning channel
dgamma_dT = (Params.delay_c1*Params.forcing_per_2xco2/Params.lambda_0);

% Two period ahead states for the value function channels, first cell is
% year 2025 state, need to get ith+1 state and ith-1 state quadratures
for i = 2:size(s,1)-2
    
    s_base_2{i-1} = [k_quad(i-1,:)' mlo_quad(i-1,:)' mup_quad(i-1,:)'...
        m_quad(i-1,:)' temp_quad(i-1,:)' to_quad(i-1,:)'...
        exp_quad(i-1,:)' var_quad(i-1,:)'];
    
    % Expected state
    s_base_exp_2{i-1} = [exp_k_2(i-1) exp_mlo_2(i-1) exp_mup_2(i-1) exp_m_2(i-1) exp_temp_2(i-1) exp_to_2(i-1) exp_exp_2(i-1) exp_var_2(i-1)];
    
    % Certainty-equivalent state, expected initial prior mean, zero
    % variance
    s_base_ce_2{i-1} = [exp_k_2(i-1) exp_mlo_2(i-1) exp_mup_2(i-1) exp_m_2(i-1) exp_temp_2(i-1) exp_to_2(i-1) Params.meanstart 0];
end

% Value function derivatives two periods ahead, column t is for time t+2
% using time t information
for t = 1:size(s,1)-3
    
    % Set correct year (current year so t = 1 implies currently at period
    % 0, year 2005)
    space.year = (t-1+1)*Params.timestep;
    
    % Regular SCC components %
    dVstar_dT_2(:,t) = smol_eval(space,[],s_base_2{1,t},t_deriv,space.dist);
    dVstar_dmu_2(:,t) = smol_eval(space,[],s_base_2{1,t},mu_deriv,space.dist);
    dVstar_dsigma_2(:,t) = smol_eval(space,[],s_base_2{1,t},var_deriv,space.dist);

end

% Value function derivatives, column t is for time t+1 using time t
% information
for t = 1:size(s,1)-2
    
    % Set correct year (current year so t = 1 implies currently at period
    % 0, year 2005)
    space.year = (t-1)*Params.timestep;
    
    % Regular SCC components, you can use these to calculate the SCC
    dVstar_dM(:,t) = smol_eval(space,[],s_base{1,t},m_deriv,space.dist);
    dVstar_dT(:,t) = smol_eval(space,[],s_base{1,t},t_deriv,space.dist);
    dVstar_dk(:,t) = smol_eval(space,[],s_base{1,t},k_deriv,space.dist);
    value(:,t) = smol_eval(space,[],s_base{1,t});
    
end

% Value function derivatives one period ahead, column t is for time t+1
for t = 1:size(s,1)-2
    
    % Set correct year (current year so t = 1 implies currently at period
    % 1, year 2015)
    space.year = (t-1)*Params.timestep;
   
    % CO2 uncertainty decomposition components %
    
    % Certainty-equivalent
    dVstar_dM_ce(:,t) = smol_eval(space,[],s_base_ce{1,t},m_deriv,space.dist);
    
    % At expected state
    dVstar_dM_exp(:,t) = smol_eval(space,[],s_base_exp{1,t},m_deriv,space.dist);
    
    % Certainty-equivalent
    dVstar_dT_ce(:,t) = smol_eval(space,[],s_base_ce{1,t},t_deriv,space.dist);
    
    % At expected state
    dVstar_dT_exp(:,t) = smol_eval(space,[],s_base_exp{1,t},t_deriv,space.dist);
    

end

% Finite difference approximation to mixed derivatives, Smolyak_Polynomial
% does not work for doing mixed derivatives
sim.diff_ep = 1e-5;
m_vec = [0 0 0 sim.diff_ep 0 0 0 0];
t_vec = [0 0 0 0 sim.diff_ep 0 0 0];
mu_vec = [0 0 0 0 0 0 sim.diff_ep 0];

% Calculate value partials for uncertainty channels one period ahead
for t = 1:size(s,1)-2
    space.year = (t-1)*Params.timestep;
    
    d3Vstar_dT3_exp(:,t) = smol_eval(space,[],s_base_exp{1,t},t_deriv_3,space.dist);
    
    d3Vstar_dM_dT2_exp(:,t) = (smol_eval(space,[],s_base_exp{1,t}+m_vec,t_deriv_2,space.dist)...
        - smol_eval(space,[],s_base_exp{1,t}-m_vec,t_deriv_2,space.dist))/(2*sim.diff_ep);
    
    d3Vstar_dmu2_dM_exp(:,t) = (smol_eval(space,[],s_base_exp{1,t}+m_vec,mu_deriv_2,space.dist)...
        - smol_eval(space,[],s_base_exp{1,t}-m_vec,mu_deriv_2,space.dist))/(2*sim.diff_ep);
    
    
    d3Vstar_dmu2_dT_exp(:,t) = (smol_eval(space,[],s_base_exp{1,t}+t_vec,mu_deriv_2,space.dist)...
        - smol_eval(space,[],s_base_exp{1,t}-t_vec,mu_deriv_2,space.dist))/(2*sim.diff_ep);
    
    d3Vstar_dT2_dmu_exp(:,t) = (smol_eval(space,[],s_base_exp{1,t}+mu_vec,t_deriv_2,space.dist)...
        - smol_eval(space,[],s_base_exp{1,t}-mu_vec,t_deriv_2,space.dist))/(2*sim.diff_ep);
    
    d3Vstar_dM_dT_dmu_exp(:,t) = (smol_eval(space,[],s_base_exp{1,t}+mu_vec+t_vec,m_deriv,space.dist)...
        - smol_eval(space,[],s_base_exp{1,t}+mu_vec-t_vec,m_deriv,space.dist)...
        - smol_eval(space,[],s_base_exp{1,t}-mu_vec+t_vec,m_deriv,space.dist)...
        + smol_eval(space,[],s_base_exp{1,t}-mu_vec-t_vec,m_deriv,space.dist))...
        /(4*sim.diff_ep^2);
end

% Active learning mean: two periods ahead

% Make new 3D array for replicated quadrature arrays
exp_next_2(1,:,:) = exp_next;
dT_de_2(1,:,:) = dT_de(1,:);

exp_next_2 = repmat(exp_next_2,Params.nqnodes.^2,1,1);
exp_next_2 = reshape(exp_next_2,size(exp_next_2,1)*size(exp_next_2,2),size(exp_next_2,3))';

dT_de_2 = repmat(dT_de_2,Params.nqnodes.^4,1,1);
dT_de_2 = reshape(dT_de_2,size(dT_de_2,1)*size(dT_de_2,2),size(dT_de_2,3))';

exp_next_2(end,:) = [];

% dmu_t+2/dT_t+1
dmu_dT = bsxfun(@rdivide,(bsxfun(@times,(bsxfun(@times,s(2:end-1,8).*dgamma_dT,H_quad)),(s(2:end-1,8).*gamma_param(2:end).^2+Params.tvar)) -...
    bsxfun(@times,(2*s(2:end-1,8).*gamma_param(2:end)*dgamma_dT),(bsxfun(@times,s(2:end-1,8).*gamma_param(2:end),H_quad)+Params.tvar*exp_next_2))),...
    (s(2:end-1,8).*gamma_param(2:end).^2 + Params.tvar).^2);


% Change in variance two periods ahead from marginal increase in emissions
% today
if LogicalOpts.fb_learn
    dsigma_de = dsigma_dt(2:end).*(dT_de(1,:));
    
    % dmu_t+2/dT_t+1 * dT_t+1/de_t
    dmu_de_2 = dmu_dT.*dT_de_2;
    
else
    dsigma_de = 0.*(dT_de(1,:));
    
    % dmu_t+2/dT_t+1 * dT_t+1/de_t
    dmu_de_2 = dmu_dT.*dT_de_2;
    dmu_de_2 = 0*dmu_de_2;
    
end

% Active learning variance: dV_t+2/dsigma_t+2 dsigma_t+2/dT_t+1 dT_t+1/de_t
channels.active_learning = (-1000*(VarsExog.betaeff(1)^Params.timestep)*...
    (weights2'*dVstar_dsigma_2).*dsigma_de(1:end-1).*VarsExog.eff_adjustn(1:end-2)'./(weights'*dVstar_dk(:,1:end-1)))';

% Active learning mean: dV_t+2/dmu_t+2 dmu_t+2/dT_t+1 dT_t+1/de_t
channels.active_learning_2 = (-1000*(VarsExog.betaeff(1)^Params.timestep)*...
    weights2'*(bsxfun(@minus,dmu_de_2(1:end-1,:)',weights2'*(dmu_de_2(1:end-1,:)')).*bsxfun(@minus,dVstar_dmu_2,weights2'*dVstar_dmu_2))...
    .*VarsExog.eff_adjustn(1:end-2)'./(weights'*dVstar_dk(:,1:end-1)))';

% Certainty-equivalent: dV_t+1/dM_t+1|ce 
ce_m = -1000*dVstar_dM_ce.*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk);
ce_t = -1000*dVstar_dT_ce.*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk).*dT_de(1,:);

% Adjustment for future uncertainty: (dV_t+1/dM_t+1|exp - dV_t+2/dM_t+1|ce)
unc_adjust_m = -1000*(dVstar_dM_exp - dVstar_dM_ce).*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk);
unc_adjust_t = -1000*(dVstar_dT_exp - dVstar_dT_ce).*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk).*dT_de(1,:);

% Precautionary abatement: 1/2 * (d3V_t+1/dM_t+1 dT_t+1^2)* var(T_t+1) 
precaution_m_1 = -1000*1/2*d3Vstar_dM_dT2_exp.*var_td(1:end-1)'.*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk);
precaution_t_1 = -1000*1/2*d3Vstar_dT3_exp.*var_td(1:end-1)'.*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk).*dT_de(1,:);

% Precautionary abatement: (d3V_t+1/dM_t+1^ dmu_t+1 dT_t+1)* cov(T_t+1,mu_t+1)
precaution_m_2 = -1000*d3Vstar_dM_dT_dmu_exp.*cov_t_e(1:end-1)'.*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk);
precaution_t_2 = -1000*d3Vstar_dT2_dmu_exp.*cov_t_e(1:end-1)'.*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk).*dT_de(1,:);

% Signal smoothing: 1/2 * (d3V_t+1/dmu_t+2^1 dM_t+1)* var(mu_t+1)
smoothing_m = -1000*1/2*d3Vstar_dmu2_dM_exp.*var_ed(1:end-1)'.*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk);
smoothing_t = -1000*1/2*d3Vstar_dmu2_dT_exp.*var_ed(1:end-1)'.*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk).*dT_de(1,:);

% Sum of Taylor expansion channels
channel_m = (ce_m + unc_adjust_m + smoothing_m + precaution_m_1 + precaution_m_2);
channel_t = (ce_t + unc_adjust_t + smoothing_t + precaution_t_1 + precaution_t_2);

% Change in temperature two periods ahead
dT_de_2 = dT_de(:,2:end);
dT_dM_de = bsxfun(@plus,bsxfun(@times, dT_dT_2(1:end-1,:)',dT_de_2(1,:)), dT_de(1,2:end)*Params.carbonxfer(1,1));

% Insurance channel two periods ahead
channels.insurance = -1000*(VarsExog.betaeff(1)^Params.timestep)*...
    (weights2'*(bsxfun(@minus,dVstar_dT_2,weights2'*dVstar_dT_2).*bsxfun(@minus,dT_dM_de,weights2'*dT_dM_de)))./(weights'*dVstar_dk(:,1:end-1));

% One period ahead partial of continuation value with respect to emissions
deriv_cont_value = dVstar_dM + dT_de.*dVstar_dT;

% Leading scaling term for EZW
scc_term_1 = (weights'*((Params.rho.*value).^(Params.risk/Params.rho))).^((Params.rho-Params.risk)/Params.risk);

% Value term in covariance for EZW
scc_term_2 = ((Params.rho.*value).^((Params.risk-Params.rho)/Params.rho));
scc_term_2_exp = (weights'*((Params.rho.*value).^((Params.risk-Params.rho)/Params.rho)));

% Covariance term for EZW
if Params.risk ~= Params.rho
    cov_val_mval = (weights'*(bsxfun(@minus,scc_term_2,scc_term_2_exp) .*...
        bsxfun(@minus,-deriv_cont_value,-weights'*deriv_cont_value)));
else
    cov_val_mval = 0;
end

% Scaling term on main SCC channel for EZW (approx = to 1)
scaling_term = (scc_term_1.*scc_term_2_exp)';

% Preference for early resolution of uncertainty term
channels.pref_temp_res = (1000.*scc_term_1.*(cov_val_mval).*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk))';

% Different channels, except temporal resolution: all in $/tC
channels.ce = scaling_term.*(ce_t' + ce_m');
channels.precaution = scaling_term.*(precaution_m_1' + precaution_m_2' + precaution_t_1' + precaution_t_2');
channels.smoothing = scaling_term.*(smoothing_m' + smoothing_t');
channels.unc_adjust = scaling_term.*(unc_adjust_m' + unc_adjust_t');
channels.total = scaling_term.*(channels.ce + channels.precaution + channels.smoothing + channels.unc_adjust);

% Main component of SCC for EZW: all in $/tC
scc_component_1  = (-1000.*...
    scaling_term'.*(weights'*deriv_cont_value).*VarsExog.eff_adjustn(1:end-1)'./(weights'*dVstar_dk))';

% Total SCC: all in $/tC
scc = scc_component_1 + channels.pref_temp_res;

% Calculate FOC errors
df_dc = VarsExog.At.^(Params.rho/(1-Params.kappa)).*VarsExog.Lt.*((qi(:,1)).^(Params.rho-1));
dk_dq1 = -Params.timestep*VarsExog.At.^(1/(1-Params.kappa)).*VarsExog.Lt./VarsExog.eff_adjustn;
dk_dq2 = -Params.timestep.*net_output./VarsExog.eff_adjustn;
dm_dmu = Params.timestep*(VarsExog.sigma.*(-1).*At.*(capital).^Params.kappa.*Lt.^(1-Params.kappa));
dmu_dq2 = (1/Params.a2).*(qi(:,2)./VarsExog.Psi).^(1/Params.a2-1)./VarsExog.Psi;
df_dm = Params.forcing_per_2xco2/log(2)./s(2:end,4);
dt_dm = Params.delay_c1.*df_dm;

error_ck = (VarsExog.betaeff(1).^Params.timestep).*(weights'*bsxfun(@times,dVstar_dk,dk_dq1(1:end-1)'));

% Consumption error
foc_error = log10(abs(-error_ck./df_dc(1:end-1)'/Params.timestep-1))';

% Abatement error
foc_error(:,2) = log10(abs(-weights'*bsxfun(@times,dVstar_dM,(dm_dmu(1:end-1).*dmu_dq2(1:end-1))')./ ...
    (weights'*bsxfun(@times,dVstar_dT,(dt_dm(1:end-1).*dm_dmu(1:end-1).*dmu_dq2(1:end-1))') +...
    weights'*bsxfun(@times,dVstar_dk,dk_dq2(1:end-1)'))-1))';

