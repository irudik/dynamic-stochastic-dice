function [VarsExog] = transitions_exogenous(Params,LogicalOpts)

% Stores the trajectories of the exogenous variables

% Time horizon of the model
horizon = Params.final_year;

% Vector of periods
t = 0:horizon/Params.timestep;

% Annual non-industrial emissions at all time nodes
for s = 0:horizon/Params.timestep
    VarsExog.Bs(s+1,:) = Params.B0*Params.gB_dec.^(s);%*Params.timestep;
end

% Labor and labor growth at each time node
VarsExog.gL(:,1) = (exp(Params.dL_dec.*(0:(horizon/Params.timestep+1)))-1)./exp(Params.dL_dec.*(0:(horizon/Params.timestep+1)));
VarsExog.Lt(:,1) = Params.L0+(Params.Linf-Params.L0)*(VarsExog.gL(1:end-1));

% Next period labor
VarsExog.Ltn(:,1) = Params.L0+(Params.Linf-Params.L0)*(VarsExog.gL(2:end));

% Technological growth rate at each time node
VarsExog.gA(:,1) = Params.gA_N*10*exp(-Params.dA*Params.timestep*(0:(horizon/Params.timestep+1)));

% Current technology
VarsExog.At(1) = Params.A0_N;
for s = 1:horizon/Params.timestep
    VarsExog.At(s+1,:) = VarsExog.At(s)./(1-VarsExog.gA(s));
end

% Next period technology
for s = 1:horizon/Params.timestep+1
    VarsExog.Atn(s,:) = VarsExog.At(s)./(1-VarsExog.gA(s));
end

% Discount factor
VarsExog.betaeff(1) = 1/(1+Params.delta);

% Exogenous decarbonization of production
VarsExog.dsigma = Params.gsigma0_dec*exp(-10*Params.dsigma*(1:horizon/Params.timestep+1)');
VarsExog.sigma(1,1) = Params.sigma0;
for s = 1:horizon/Params.timestep
    VarsExog.sigma(s+1,:) = VarsExog.sigma(s)./(1-VarsExog.dsigma(s));
end

% Emission intensity of production, adjust technology so no longer labor-augmenting
VarsExog.emint = VarsExog.sigma.*VarsExog.At.^Params.kappa.*VarsExog.Lt;

% Abatement cost function coefficient at all time nodes
VarsExog.Psi = (Params.a0*VarsExog.sigma./Params.a2).*((Params.a1 - 1+ exp(-Params.gPsi_dec*(0:horizon/Params.timestep)' ))/Params.a1);

% External forcing
VarsExog.EF = (Params.EF0 + 0.1*(Params.EF1 - Params.EF0)*min(t,10))';

% External forcing in next year
VarsExog.EF_next = (Params.EF0 + 0.1*(Params.EF1 - Params.EF0)*min(t+1,10))';

% Adjustments to put capital and consumption into effective terms
% Dollars per transformed-technology and labor
% This allows for a better quality approximation, and quicker solve
VarsExog.eff_adjust = VarsExog.At.^(1/(1-Params.kappa)).*VarsExog.Lt;
VarsExog.eff_adjustn = VarsExog.Atn.^(1/(1-Params.kappa)).*VarsExog.Ltn;

end

