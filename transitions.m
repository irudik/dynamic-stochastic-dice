function [snext, maxo, mino] = transitions(si,~,consumption,abatecost,...
    vars,Params,LogicalOpts,space,feedback,capital1)

%   transitions: Returns next periods state conditional on the control that
%   maximized the right hand side of the Bellman. Returns the maximum and
%   minimum evaluation points for surface temperature and the mean belief
%   due to the quadrature for the feedback parameter.

% Next state vector
snext = zeros(1,length(si));

% Current abatement
mu = (abatecost/vars.Psii)^(1/Params.a2);

% State dependent exogenous variables
capital = si(1,1)*vars.eff_adjusti;
output = vars.Ati*(capital)^Params.kappa*vars.Lti^(1-Params.kappa);
net_output = vars.Ati*(capital)^Params.kappa*vars.Lti^(1-Params.kappa)./(1+Params.b2*si(1,5).^Params.b3);

% Initialize noise shocks
noise_a = 0;
noise_o = 0;


% Capital
snext(:,1) = ((capital)*(1-Params.deltak)^Params.timestep + ...
    Params.timestep*(net_output*(1-abatecost) - vars.eff_adjusti*consumption))/...
    vars.eff_adjustni;

% CO2 - Lower Ocean
snext(:,2) = Params.carbonxfer(3,3)*si(1,2) + Params.carbonxfer(2,3)*si(1,3);

% CO2 - Upper Ocean
snext(:,3) = Params.carbonxfer(2,2)*si(1,3) + Params.carbonxfer(1,2)*si(1,4) + ...
    Params.carbonxfer(3,2)*si(1,2);

% CO2 - Atmosphere
snext(:,4) = Params.timestep*(vars.sigmai*(1-mu)*output + vars.Bi) +...
    Params.carbonxfer(1,1)*si(1,4) + Params.carbonxfer(2,1)*si(1,3);

forcing = Params.forcing_per_2xco2*(log((snext(:,4))/Params.Mpre)/log(2))+vars.EFi_next;

% Surface temperature at the real cs
snext(1,5) = si(:,5) + ...
    Params.delay_c1.*(forcing - (Params.forcing_per_2xco2*(1-Params.fb_cert)/(Params.lambda_0))*si(:,5) + ...
    Params.delay_c3*(si(:,6)-si(:,5))) + noise_a;

% Surface temperature at all quadrature nodes
temp_next = (si(:,5) + ...
    Params.delay_c1.*(forcing - (Params.forcing_per_2xco2*(1-feedback)/(Params.lambda_0))*si(:,5) + ...
    Params.delay_c3*(si(:,6)-si(:,5))))';% + noise_a;

temp_comb = combvec(temp_next,Params.te')';

% Get surface temperature quadrature nodes
temp_nodes = temp_comb(:,1) + temp_comb(:,2);

% Maximum and minimum temperature quadrature nodes
maxo = max(temp_nodes);
mino = min(temp_nodes);

% Ocean temperature
snext(:,6) = (1-Params.delay_c4)*si(:,6) + Params.delay_c4*si(:,5);

if LogicalOpts.fb_learn % Learning
    
    % H and gamma that correspond to the paper
    H = snext(:,5) - (si(1,5) + Params.delay_c1.*(forcing ...
        - (Params.forcing_per_2xco2./(Params.lambda_0))*si(1,5) + Params.delay_c3*(si(1,6)-si(1,5))));
    
    gamma = (Params.delay_c1*Params.forcing_per_2xco2*si(1,5)/Params.lambda_0);
    
    % Next period's feedback mean
    snext(:,7) = (si(1,8)*(gamma)*(H) + Params.tvar*si(1,7))/...
        (Params.tvar+si(1,8)*(gamma).^2);
    
    % Next period's feedback variance
    snext(:,8) = Params.tvar*si(1,8)/(Params.tvar+si(1,8)*(gamma).^2);
    
    H = temp_nodes - (si(1,5) + Params.delay_c1.*(forcing ...
        - (Params.forcing_per_2xco2./(Params.lambda_0))*si(1,5) + Params.delay_c3*(si(1,6)-si(1,5))));
    
    % Get mean belief quadrature nodes
    exp_nodes = (si(1,8)*(gamma)*(H) + Params.tvar*si(1,7))/...
        (Params.tvar+si(1,8)*(gamma).^2);
    
    % Maximum and minimum mean belief quadrature nodes
    maxo(2) = max(exp_nodes);
    mino(2) = min(exp_nodes);
    
else % No learning
    
    % Next period's feedback mean
    snext(:,7) = si(1,7);
    
    % Next period's feedback variance
    snext(:,8) = si(1,8);
    
    maxo(2) = snext(:,7);
    mino(2) = snext(:,7);
end

end