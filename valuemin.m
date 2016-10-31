function [f,g,H] = valuemin(q,si,vars,Params,LogicalOpts,space,...
    space_terminal,feedback,weights,gradient)

%   valuemin: Returns the right hand side of the Bellman conditional on the
%   effective consumption and abatement cost values selected by the solver.

if size(q,1) == 2
    q = q'; %convert 2x1 q to 1x2
end

% Next period's state vector
g = zeros(length(feedback)*length(Params.te),length(si));

% Abatement rate
mu = (q(1,2)/vars.Psii)^(1/Params.a2);

% Next period's atmospheric CO2 stock
g(:,4) = Params.timestep*(vars.sigmai*(1-mu)*vars.output + vars.Bi) +...
    Params.carbonxfer(1,1)*si(1,4) + Params.carbonxfer(2,1)*si(1,3);

% Next period's upper ocean CO2 stock
g(:,3) = Params.carbonxfer(2,2)*si(1,3) + Params.carbonxfer(1,2)*si(1,4) + ...
    Params.carbonxfer(3,2)*si(1,2);

% Next period's lower ocean CO2 stock
g(:,2) = Params.carbonxfer(3,3)*si(1,2) + Params.carbonxfer(2,3)*si(1,3);

% Next period's forcing
forcing = Params.forcing_per_2xco2*(log(g(1,4)/Params.Mpre)/log(2))+vars.EFi_next;

% Next period's surface temperature, feedback is the quadrature for the
% feedback parameter
temp_next = (si(1,5) + Params.delay_c1.*(forcing ...
    - (Params.forcing_per_2xco2*((1-feedback)./(Params.lambda_0)))*si(1,5) + Params.delay_c3*(si(1,6)-si(1,5))))';

% Get all possible combinations of temperature quadrature (without the
% shock) and the temperature shock quadrature (Total of Params.nqnodes^2
% quadrature nodes for the 2D integration)
temp_comb = combvec(temp_next,Params.te')';

% Temperature + temperature shock
g(:,5) = temp_comb(:,1) + temp_comb(:,2);

% Next period's ocean temperature
g(:,6) = ((1-Params.delay_c4)*si(1,6) + Params.delay_c4*si(1,5));

% Next period's capital stock in effective terms
g(:,1) = ((vars.capital)*(1-Params.deltak)^Params.timestep + ...
    Params.timestep*(vars.net_output*(1-q(1,2)) - vars.eff_adjusti*q(1,1)))/...
    vars.eff_adjustni;

if LogicalOpts.fb_learn
    
    % H and gamma from the paper
    H = g(:,5) - (si(1,5) + Params.delay_c1.*(forcing ...
    - (Params.forcing_per_2xco2./(Params.lambda_0))*si(1,5) + Params.delay_c3*(si(1,6)-si(1,5))));

    gamma = (Params.delay_c1*Params.forcing_per_2xco2*si(1,5)/Params.lambda_0);
    
    % Next period's feedback mean
    g(:,7) = (si(1,8)*(gamma)*(H) + Params.tvar*si(1,7))/...
        (Params.tvar+si(1,8)*(gamma).^2);

    % Next period's feedback variance
    g(:,8) = Params.tvar*si(1,8)/((Params.tvar+si(1,8)*(gamma).^2));
    
else
    
    % Next period's feedback mean
    g(:,7) = si(1,7);

    % Next period's feedback variance
    g(:,8) = si(1,8);
    
end

% Evaluate value function at future state quadrature
V = smol_eval(space,space_terminal,g);

% If negative continuation value, it should be positive (not working yet)
if Params.rho > 0 && any(V < 0)
     V = abs(V);
end

% Continuation value
V = (V*Params.rho).^(Params.risk).^(1/Params.rho);
Ev = weights'*V(:,1);

% Utility function
f = vars.Ati^(Params.rho/(1-Params.kappa))*vars.Lti*((q(1,1)).^Params.rho)/Params.rho;

% The right hand side of the Bellman equation
cont_value = (vars.betaeffi^Params.timestep)*(Ev.^(Params.rho/Params.risk))/Params.rho;

f = -Params.timestep*f - cont_value;

% avoid numerical infinity
if f>realmax
    f=realmax;
end

% Return gradient and hessian if asked
if nargout  > 1
    g = gradient(q);
end

end
