clear;
close all;
clc;

set(groot, 'defaultTextInterPreter', 'tex')
set(groot, 'defaultAxesTickLabelInterPreter', 'tex')
set(groot, 'defaultFigureColormap', turbo(256));
set(groot, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');
set(groot, 'defaultSurfaceEdgeAlpha', 0.3);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultFigureColor', [1; 1; 1]);
set(groot, 'defaultAxesColor', 'none');
set(groot,'DefaultAxesFontSize', 15);

%% PATH GENERATION

if isunix()
    CEApath = "CEA/";
else
    CEApath = "CEA\";
end

addpath(genpath(CEApath));

%% ENGINE DATA

engine.t_b = 50; % Burning time
engine.P_c = 50; % Combustion chamber pressure
utils.g = 9.81; 
utils.R = 8.3145; 
utils.Rair = 287.2; % Specific universal constant
engine.eps = 20; % Expansion ratio
engine.of = 2.24; % OF ratio
engine.dmp = 100; % Propellant mass flow 
ox.dm = engine.dmp * (engine.of/(engine.of+1)); % oxidier mass flow
f.dm = 1/(engine.of+1) * engine.dmp; % fuel mass flow
ox.m = ox.dm * engine.t_b; % oxidier mass
f.m = f.dm * engine.t_b; % fuel mass
f.rho = 810; % Fuel density
ox.rho = 1140; % Oxidier density
ox.V = ox.m/ox.rho; % Oxidier volume
f.V = f.m/f.rho; % Fuel volume

%% CEA
% see how it works

ceaOutout = CEA('problem','rocket','frozen','nfz',1,'p,pa',engine.P_c,'o/f', ...
engine.of, ...
    'reactants','fuel','RP-1(L)', 'C',1.,'H',1.9423,...
    'oxid','O2(L)','end');
%% PERFORMANCES

% Data input from CEA simulation
perf.M = ceaOutout.output.froz.mach(end); % Mach in CC
perf.Ve_i = ceaOutout.output.froz.sonvel(end) * perf.M; % Flow velocity
perf.isp = ceaOutout.output.froz.isp(end); % Specific impulse
perf.ispVac = ceaOutout.output.froz.isp_vac(end); % Specific impulse in vacuum
gamma = ceaOutout.output.froz.gamma(2); % specific heat ratio
perf.gamma = gamma;
perf.cS = ceaOutout.output.froz.cstar(1); % Characteristic Velocity (c*)
perf.cT_i = ceaOutout.output.froz.cf(end); % Thrust coefficient (cT)
perf.Pe = ceaOutout.output.froz.pressure(end); % Exhaust Pressure (P_e)
perf.Tc = ceaOutout.output.froz.temperature(1);
perf.Tt = ceaOutout.output.froz.temperature(2);
perf.Pt = ceaOutout.output.froz.pressure(2);

% Application of 2D loss coefficient
nozzle.alphaDiv = deg2rad(15);
nozzle.lambda = (1+cos(nozzle.alphaDiv))/2;

perf.Ve = perf.Ve_i * nozzle.lambda;

%% AREAS AND NOZZLE

nozzle.At = perf.cS * engine.dmp / engine.P_c;
nozzle.Dt = sqrt(4*nozzle.At/pi);
nozzle.Ae = nozzle.At * engine.eps;
nozzle.De = sqrt(4*nozzle.Ae/pi);

%% THRUST

% Time vector
time = 0:0.1:engine.t_b;

% Atmospheric Pressure 
altRange = linspace(0,30000,length(time));
[~, ~, Patm, ~] = atmosisa(altRange,extended=true);

% Thrust 

perf.T = engine.dmp * perf.Ve + (perf.Pe - Patm) * nozzle.Ae;


%% PERFORMANCES 2

perf.cT = engine.dmp * perf.Ve / (engine.P_c  * nozzle.At);
perf.itot = trapz(time,perf.T);
perf.ispVec = perf.T / engine.dmp / 9.81;
engine.rhoAve = (ox.m + f.m) / (ox.V + f.V);
perf.ivVec = perf.ispVec * engine.rhoAve;
perf.iv = perf.isp * engine.rhoAve;


%% AREAS AND NOZZLE

nozzle.Dt = sqrt(4*nozzle.At/pi);

machFun = @(m) 1./m .* (2/(gamma+1) .* (1 + (gamma-1)/2 .* m.^2)) .^ ((gamma+1)/(2*(gamma-1)));

nozzle.Ac = nozzle.At * machFun(0.2);
nozzle.Dc = sqrt(4*nozzle.Ac/pi);

nozzle.alphaDiv = deg2rad(15);
nozzle.alphaConv = deg2rad(45);
nozzle.lConv = (nozzle.Dc-nozzle.Dt) / (2 * tan(nozzle.alphaConv));
nozzle.Ae = engine.eps * nozzle.At;
nozzle.De = sqrt(4*nozzle.Ae/pi);
nozzle.lDiv = (nozzle.De-nozzle.Dt) / (2 * tan(nozzle.alphaDiv));
nozzle.r = [linspace(nozzle.Dc/2,nozzle.Dt/2,500), linspace(nozzle.Dt/2,nozzle.De/2,2000)];
nozzle.areas = pi * nozzle.r.^2;
nozzle.longitudinalCoordinate = [linspace(0,nozzle.lConv,500),linspace(nozzle.lConv,nozzle.lDiv,2000)];

%% ISOENTROPIC EXPANSION IN THE NOZZLE

guess = [linspace(0.2,1,500), linspace(1,5,2000)];

solveMachFun = @(m) machFun(m) - nozzle.areas / nozzle.At;

engine.mach = fsolve(solveMachFun,guess);

nozzle.areaRatios = machFun(engine.mach);


engine.temp = (1 + (gamma-1)/2 * engine.mach.^2).^-1 .* ceaOutout.output.froz.temperature(1); 


mu = 1.0698*1e-4;
Pr = 0.5512;
k_c = 15.3748*1e-3/1e-2;
engine.bartz.recoveryFactor = Pr^(1/3); %recovery factor 

Twi = (-150) +273.15;

Tcc_adiabatic = engine.temp(1) * (1 + engine.bartz.recoveryFactor * (ceaOutout.output.froz.gamma(1)-1)/2 * 0.2);

sWall = 1.5e-3; %wall thickness
K_copper= 390; %[W/m*K]
engine.bartz.cpProducts = ceaOutout.output.froz.cp(2) * 1e3;


engine.bartz.sigma = @(Tw) 1./(( 1/2 * Tw/engine.temp(1) .* (1+(gamma-1)/2 .* engine.mach.^2) + 1/2 ).^ (0.8-0.6/5) .* (1+(gamma-1)/2 * engine.mach.^2).^(0.6/5));

h = @(Tw) ( (0.026/nozzle.Dt^0.2) * ((mu^0.2*engine.bartz.cpProducts)/Pr^0.6) * (engine.P_c/ceaOutout.output.froz.cstar(end))^0.8 ) .* (1./nozzle.areaRatios).^0.9 .* engine.bartz.sigma(Tw);

qWall = @(Tw) K_copper / sWall * (Tw-Twi);


eq = @ (Tw) h(Tw) .* (Tcc_adiabatic-Tw) - qWall(Tw);

guess = 300 * ones(1,length(nozzle.areaRatios));

Twall = fsolve(eq,guess);


%% PLOTS and EXPORT

% Regenerative cooling
figure();
plot(nozzle.longitudinalCoordinate,Twall,'LineWidth',5)
ylabel("\textbf{Inside wall Temperature [K]}",Interpreter="latex",FontSize=40);
xlabel("\textbf{Longitudinal coordinate [m]}",Interpreter="latex",FontSize=40);
xline(2.5133,'k--',Label="End of regenrative cooling",Interpreter="latex",FontSize=20,LineWidth=3)

yyaxis right
plot(nozzle.longitudinalCoordinate,nozzle.r,'LineWidth',5)
ylabel("\textbf{Nozzle radius [m]}",Interpreter="latex",FontSize=40);

legend("Twall",FontSize=35)
axis padded
grid minor
title("Regenerative cooling");

% Thrust and Pressure over time
figure();
plot(time,perf.T,'LineWidth',3)
xlabel("\textbf{Time [s]}",Interpreter="latex",FontSize=20);
ylabel("\textbf{Thrust [N]}",Interpreter="latex",FontSize=20);
hold on

yyaxis right
plot(time,Patm,'LineWidth',3);
ylabel("\textbf{Ambient Pressure [Pa]}",Interpreter="latex",FontSize=20);

legend("Thrust")
axis padded
grid minor
title("Thrust and Pressure over time");


% Plot
figure();
plot(time,perf.ispVec,'LineWidth',3)
yline(perf.isp,'k--','LineWidth',2);
xlabel("\textbf{Time [s]}",Interpreter="latex",FontSize=20);
ylabel("\textbf{Specific Impulse [s]}",Interpreter="latex",FontSize=20);

yyaxis right
plot(time,Patm,'LineWidth',3);
ylabel("\textbf{Ambient Pressure [Pa]}",Interpreter="latex",FontSize=20);

legend("Specific Impulse","OE Specific Impulse","AutoUpdate","off",'Location','best');
axis padded
grid minor
title("Specific Impulse and Pressure over time");

% Mach number
figure();
plot(nozzle.longitudinalCoordinate,engine.mach,'LineWidth',2)
ylabel("\textbf{Mach [-]}",Interpreter="latex",FontSize=20);
xlabel("\textbf{Longitudinal coordinate [m]}",Interpreter="latex",FontSize=20);

yyaxis right
plot(nozzle.longitudinalCoordinate,nozzle.r,'LineWidth',2)
ylabel("\textbf{Nozzle radius [m]}",Interpreter="latex",FontSize=20);

legend("Mach")
axis padded
grid minor
title("Mach number");


