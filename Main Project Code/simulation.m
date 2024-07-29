clear all
close all
%clc

%% Design parameters

eFeedLines = 2e-6; % Feed line roughness [m]
eCoolingJacket = 50e-6; % Cooling jacket roughness [m]

TotalVolume = 0.1;%1.2565;
OxDensity = 980;
FuDensity = 830;
ocup = 0.458;
of = 2.4*FuDensity/OxDensity;
oxVolume = TotalVolume*(of)/(of+1);
fuVolume = TotalVolume*(1)/(of+1);
mOx = oxVolume*ocup*OxDensity;
mFu = fuVolume*ocup*FuDensity;
% Oxidizer Tank

oxTank.volume = oxVolume; % Volume of the oxidizer tank [m^3]
oxTank.initialPress = 5.95e+06; % Initial pressure of the oxidizer tank [Pa]
oxTank.initialTemp = 120; % Initial temperature of oxidizer and its pressurant [K]
oxTank.feedLength = 0.83; % Length of the oxidizer feed line [m]
oxTank.feedDiameter = 0.00625; % Diameter of the oxidizer feed line [m]
oxTank.initialMass = mOx; % Initial mass of oxidizer [kg]

% Fuel Tank

fuTank.volume = fuVolume; % Volume of the fuel tank [m^3]
fuTank.initialPress = 6.1e+06; % Initial pressure of the fuel tank [Pa]
fuTank.initialTemp = 293; % Initial temperature of fuel and its pressurant [K]
fuTank.feedLength = 0.38; % Length of the oxidizer feed line [m]
fuTank.feedDiameter = 0.00625; % Diameter of the fuel feed line [m]
fuTank.initialMass = mFu; % Initial mass of fuel [kg]

% Propellant couple (RP-1\LOX)
ox.density = 980; % Density of the oxidant [kg/m^3]
ox.viscosity = 1.0280e-04; % Dynamic viscosity of the oxidant [kg/m^3]
fu.density = 830; % Density of the fuel [kg/m^3]
fu.viscosity = 7.5000e-04; % Dynamic viscosity of the fuel [kg/m^3]
fu.Cp = 1.8841e+03; % Cp of the fuel [J/kg K]
fu.prandtl = 12.2876; % Prandtl of the fuel
fu.k = 0.115; % Thermal conductivity of the fuel [W/m K]
pressurant.gamma = 1.66; % Gamma of the pressurant

% Injector Plate
injector.holeSizeOx = 0.0012; % Size of the oxidizer injector holes [m]
injector.holeSizeFu = 0.0010; % Size of the fuel injector holes [m]
injector.holeNumberOx = 8; % Number of oxidizer injector holes
injector.holeNumberFu = 4; % Number of fuel injector holes
injector.CdOx = 0.65; % Discharge coeficient for the oxidizer injector
injector.CdFu = 0.65; % Discharge coeficient for the fuel injector

% Engine
engine.combustionChamberLength = 0.127; % Length of the combustion chamber [m]
engine.combustionChamberDiameter = 0.035; % Diameter of the combustion chamber [m]
engine.throatDiameter = 0.0117; % Diameter of the throat [m]
engine.nozzleLength = 1; % Length of the diverging section [%]
engine.expansionRatio = 65; % Expansion ratio
engine.sigma_e = 4.13; % Nozzle exit angle [degrees]
engine.sigma_n = 40.18; % Nozzle initial angle [degrees]
engine.throatArea = pi*engine.throatDiameter^2/4; % Engine throat area [m^2]
engine.nPoints = 100; % Number of sections to divide the nozzle
engine.coolingWallThickness = 0.001; % Thickness of the chamber wall [m]
engine.materialConductivity = 41.9; % Condutivity of the engine material [W/mK]
engine.channelThickness = 0.005;
engine.channelThicknessT = 0.001; % Thickness of the throat cooling jacket [m]
engine.channelThicknessN = 0.003; % Thickness of the nozzle cooling jacket [m]
engine.channelThicknessC = 0.002; % Thickness of the chamber cooling jacket [m]
engine.coolingChannelNumberT = 9; % Number of throat cooling channels
engine.coolingChannelNumberN = 5; % Number of nozzle cooling channels
engine.coolingChannelNumberC = 6; % Number of chamber cooling channels
engine.throatBegin = 0.13;
engine.throatEnd = 0.145;
engine.coolingDistance = 0.25; % How far down into the engine the cooling jacket starts [m]

[engine.xPoints, engine.yPoints, engine.totalLength, engine.length, engine.alpha] = RAO_profile(engine);

%% Simulation

deltaTime = 10; % Time step [s]

% State variables

massOx = [oxTank.initialMass]; % Remaining mass of oxidizer [kg]
massFu = [fuTank.initialMass]; % Remaining mass of fuel [kg]
massFlowRateOx = [0.2313]; % Mass flow rate of oxidizer [kg/s]
massFlowRateFu = [0.0875]; % Mass flow rate of fuel [kg/s]
pressOx = [oxTank.initialPress]; % Pressure of the oxidizer tank [Pa]
pressFu = [fuTank.initialPress]; % Pressure of the fuel tank [Pa]
chamberPress = [(pressOx + pressFu)*0.85/2 - 5*10^5]; % Combustion chamber pressure [Pa]
time = [0]; % Time [s]
thrust = [0]; % Thrust [N]
throatReynolds = [3.8739e+05];
chamberGamma = [1.13];
throatGamma = [1.21];

tempFu = [fuTank.initialTemp + 300];
tempWall = ones([1,engine.nPoints-1])*fuTank.initialTemp + 300;
hgas = zeros([1,engine.nPoints-1]);

chamberTemp = [0];
%%
while massOx(end) > 0 && massFu(end) > 0
    % Elapsed time
    time(end+1) = time(end) + deltaTime;
    
    % Calculate current pressure of the oxidizer
    pressOx(end+1) = calculateIsoentropicPressure(oxTank.initialMass, oxTank.initialPress, massOx(end), ox.density, pressurant.gamma, oxTank.volume);
    
    % Calculate current pressure of the fuel
    pressFu(end+1) = calculateIsoentropicPressure(fuTank.initialMass, fuTank.initialPress, massFu(end), fu.density, pressurant.gamma, fuTank.volume);
    
    % Calculate mass flow rate of fuel and oxidizer and combustion chamber pressure
    [massFlowRateOx(end+1), massFlowRateFu(end+1), chamberPress(end+1), chamberTemp(end+1),...
        chamberGamma(end+1), throatGamma(end+1), tempFu(end+1), tempWall(end+1,:), throatReynolds(end+1), hgas(end+1,:)] = calculateChamberPress(...
        pressOx(end), pressFu(end), massFlowRateOx(end), massFlowRateFu(end), ox, fu, oxTank, fuTank,...
        eFeedLines, eCoolingJacket, injector, engine, tempWall(end,:), tempFu(end),chamberPress(end), throatReynolds(end), chamberGamma(end), throatGamma(end));
    
    % Calculate produced thrust
    thrust(end+1) = calculateThrust(chamberPress(end), engine.throatArea, throatGamma(end), engine.expansionRatio, engine.sigma_e, engine.alpha);
    
    % Calculate remaining mass of oxidizer
    massOx(end+1) = massOx(end) - massFlowRateOx(end)*deltaTime;
    
    % Calculate remaining mass of fuel
    massFu(end+1) = massFu(end) - massFlowRateFu(end)*deltaTime;
    
    % Verify fuel temperature
    %tempFu(end+1) = calculateFuelTemperature();
end

% Add some analysis on Mach number, Specific Impulse could be interesting
%% Graphs

createPlots;
