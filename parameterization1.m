%   parameterization1: Initializes the parameters of the model.

%% Time state variable

% Parameters for model period length
Params.final_year = final_year - 2005;
space.final_year = Params.final_year;
Params.terminal_value = terminal_value;
Params.timestep = 10;

% Switch for indicating we're in the final year of the horizon
Params.terminal_switch = 1;

%% Base case model

% Climate sensitivity
Params.cscert = 3;  % DICE-2007 uses 3 for climate sensitivity
Params.lambda_0 = 1.2; % Roe and Baker zero feedback climate sensitivity

% Parameters of damage function:
Params.b1 = 0; % marginal damage intercept
Params.b2 = 0.0028388; % coefficient of quadratic damage
Params.b3 = 2; % damage exponent

% Emissions from land use change and forestry
Params.B0 = 1.1; % non-industrial emissions at t=0 (GtC/yr)
Params.gB = -0.01; % growth rate of non-industrial emissions
Params.gB_dec = 0.9; % growth rate of non-industrial emissions

% Parameters of abatement function
Params.a0 = 1.17; % cost of backstop in 2005 ($1000 per ton C)
Params.a1 = 2; % ratio of initial over final backstop cost

% Parameters for the abatement cost function
Params.a2 = 2.8; % cost exponent
Params.gPsi_dec = 0.05; % discrete time growth rate of abatement cost

% Pure rate of time preference (annual)
Params.delta = 0.015;

% Parameters of population dynamics
Params.L0 = 6514; % population in 2005 (millions)
Params.dL_dec = 0.35; % growth rate of population per year (rate of convergence to asymptotic population)
Params.Linf = 8600; % asymptotic population (millions) (read inf as infinity)

% Parameters for production
Params.kappa=0.3; % capital elasticity in production

% _N=Nordhaus, who does not use A as labor-augmenting technological
% progress but instead has A L^(1-kappa) K^kappa
Params.A0_N = 0.02722; % initial production tech
Params.gA_N = 0.0092; % annual growth rate of production tech
Params.dA = 0.001; % annual rate of decline of production tech grwoth rate (*100 gives %) % Note 10-02: changed from 0.01 to 0.001: Nordhaus uses delta=0.001 as decline of tech per decade but then multiplies by 10 when calculating

% Renormalized values to labor-augmenting technological progress
% (AL)^(1-kappa) K^kappa
Params.A0 = Params.A0_N^(1/(1-Params.kappa));
Params.gA0 = Params.gA_N/(1-Params.kappa);

% Emissions intensity of production
Params.sigma0 = 0.13418; % CO2 equivalent emission-GDP ratio
Params.gsigma0_dec = -0.0730; % initial growth rate of decarbonization
Params.dsigma = 0.003; % rate of decline of decarbonization

% Annual depreciation rate of capital stock
Params.deltak = 0.1; 

% External forcing
Params.EF0 = -0.06; % exogenous forcing in 2000 (non-CO2 GHGs)
Params.EF1 = 0.30; % exogenous forcing in 2100 (non-CO2 GHGs)

%% Temperature function parameters

Params.Mpre = 596.4; % preindustrial CO2 concentration in Gt C (596.4 GtC is 276 ppm CO2)

% unit conversion parameter
Params.gtc_per_ppm = 2.16;

% forcing from 2xCO2 (DICE-2007)
Params.forcing_per_2xco2 = 3.8; % W/m2

% Temperature delay coefficients
Params.delay_c1 = 0.22; % fraction of total forcing that increments temperature
Params.delay_c3 = 0.3; % converts ocean temperature to atmospheric forcing
Params.delay_c4 = 0.05; % fraction of temperature gradient that increments ocean temperature

% DICE-2007 carbon transition model
Params.carbonxfer_dice = [81.0712 18.9288 0; 9.7213  85.2787 5; 0 0.3119 99.6881]; % from DICE-2007; element i,j is fraction from carbon pool i going to carbon pool j, meaning rows add to 100
Params.carbonxfer = Params.carbonxfer_dice./100;  % Convert to fraction

% Temperature shock parameters
Params.tvar = 0.11; % Variance of weather shocks decadally Kelly and Tan 2015, Kelly and Kolstad 1999
Params.texp = 0;

% Number of nodes for double integration
space.nqnodes = Params.nqnodes^2;

% Weather shock quadrature
try
    [Params.te,Params.tw] = qnwnorm(Params.nqnodes,Params.texp,Params.tvar);
catch
end

% Mean and variance of distribution over the feedback parameter
Params.fb_cert = 0.6;
Params.fb_var = 0.13^2;

% Turn off noise if doing a deterministic run
if ~LogicalOpts.fb_unc
    Params.te = 0*Params.te;
end

% Starting values for time paths
Params.kstart = 137/(Params.A0*Params.L0);% effective capital
Params.Mustart = 1255;  % upper ocean CO2 stock
Params.Mlstart = 18365;  % lower ocean CO2 stock
Params.Mstart = 808.9;  % atm CO2 stock
Params.Tstart = 0.7307; % surface temperature
Params.Tostart = 0.0068; % ocean temperature
Params.meanstart = Params.fb_cert; % feedback belief
Params.varstart = Params.fb_var; % feedback variance