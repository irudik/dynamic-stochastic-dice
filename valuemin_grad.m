function f = valuemin_grad(q,si,vars,Params,LogicalOpts,space,space_terminal,feedback,weights)

% valuemin_grad: The gradient of the value function.

if size(q,1) == 2
    q = q'; %convert 2x1 q to 1x2
end

% Next period's state vector, expands by 3 to capture more histories
g = zeros(length(feedback)*length(Params.te),length(si));

% Implied abatement level
mu = (q(1,2)/vars.Psii)^(1/Params.a2);

% dabatement rate/dabatement cost
dmu_dq2 = (1/Params.a2)*(q(1,2)/vars.Psii)^(1/Params.a2-1)/vars.Psii;

% Next period's atmospheric CO2 stock
g(:,4) = Params.timestep*(vars.sigmai.*(1-mu).*vars.output + vars.Bi) +...
    Params.carbonxfer(1,1)*si(1,4) + Params.carbonxfer(2,1)*si(1,3);

% dM/dabatement rate
dm_dmu = Params.timestep*(vars.sigmai.*(-1).*vars.output);

% Next period's upper ocean CO2 stock
g(:,3) = Params.carbonxfer(2,2)*si(1,3) + Params.carbonxfer(1,2)*si(1,4) + ...
    Params.carbonxfer(3,2)*si(1,2);

% Next period's lower ocean CO2 stock
g(:,2) = Params.carbonxfer(3,3)*si(1,2) + Params.carbonxfer(2,3)*si(1,3);

% Next period's forcing
forcing = Params.forcing_per_2xco2*(log(g(1,4)/Params.Mpre)/log(2))+vars.EFi_next;

% dForcing/dM
df_dm = Params.forcing_per_2xco2/log(2)/g(1,4);

% Next period's temperature, feedback is the quadrature for the feedback
% parameter
temp_next = (si(1,5) + Params.delay_c1.*(forcing ...
    - (Params.forcing_per_2xco2*((1-feedback)./(Params.lambda_0)))*si(1,5) + Params.delay_c3*(si(1,6)-si(1,5))))';
temp_comb = combvec(temp_next,Params.te')';
g(:,5) = temp_comb(:,1) + temp_comb(:,2);

% dT/dM
dt_dm = Params.delay_c1.*df_dm;

% Next period's ocean temperature plus ocean shock
g(:,6) = ((1-Params.delay_c4)*si(1,6) + Params.delay_c4*si(1,5));

% get next period's capital stock, deterministic fct of today's states +
% controls
g(:,1) = ((vars.capital)*(1-Params.deltak)^Params.timestep + ...
    Params.timestep*(vars.net_output*(1-q(1,2)) - vars.Ati^(1/(1-Params.kappa))*vars.Lti*q(1,1)))/...
    vars.eff_adjustni;



if LogicalOpts.fb_learn
    
    % H and gamma from the paper
    H = g(:,5) - (si(1,5) + Params.delay_c1.*(forcing ...
    - (Params.forcing_per_2xco2./(Params.lambda_0))*si(1,5) + Params.delay_c3*(si(1,6)-si(1,5))));
    gamma = (Params.delay_c1*Params.forcing_per_2xco2*si(1,5)/Params.lambda_0);
    
    % Next period's feedback mean
    g(:,7) = (si(1,8)*(gamma)*(H) + Params.tvar*si(1,7))/...
        (Params.tvar+si(1,8)*(gamma).^2);
    g(:,8) = Params.tvar*si(1,8)/(Params.tvar+si(1,8)*(gamma).^2);
    
else
    
    % Next period's feedback mean
    g(:,7) = si(1,7);

    % Next period's feedback variance
    g(:,8) = si(1,8);
end

% Evaluate value function at future state
dEv_dk = smol_eval(space,space_terminal,g,[1 0 0 0 0 0 0 0],space.dist);
dEv_dk = weights'*dEv_dk(:,1);
dEv_dm = smol_eval(space,space_terminal,g,[0 0 0 1 0 0 0 0],space.dist);
dEv_dm = weights'*dEv_dm(:,1);
dEv_dt = smol_eval(space,space_terminal,g,[0 0 0 0 1 0 0 0],space.dist);
dEv_dt = weights'*dEv_dt(:,1);

% Derivative of capital with respect to consumption and abatement cost
dk_dq1 = -Params.timestep*vars.Ati^(1/(1-Params.kappa))*vars.Lti/ vars.eff_adjustni;
dk_dq2 = - Params.timestep*vars.net_output/ vars.eff_adjustni;

% Marginal utility
df_dc = vars.Ati^(Params.rho/(1-Params.kappa))*vars.Lti*((q(1,1)).^(Params.rho-1));

% Consumption component of gradient
f(1) = -Params.timestep*df_dc - (vars.betaeffi^Params.timestep)*(dEv_dk)*(dk_dq1);

% Abatement cost component of gradient
f(2) = -(vars.betaeffi^Params.timestep)*(dEv_dm*dm_dmu*dmu_dq2 + ...
    dEv_dt*dt_dm*dm_dmu*dmu_dq2 + dEv_dk*dk_dq2);

end
    