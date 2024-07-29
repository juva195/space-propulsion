clear 
close all
clc

addpath('nozzle\');
%% INPUT DATA
% Constants
R = 8.3143; % universal gas costant (J/(K*mol))
g0 = 9.81; % gravity acceleration (m/s^2)
Mm1 = 21.9e-3; % Molar mass (Kg/mol)
Mm2 = 23.3e-3; % Molar mass (Kg/mol)
gamma = 1.24; % Specific heat ratio of the mixture

% Given data
T_in = 1000; % Initial Thrust (N)
P_cc_in = 50; %Combustion chamber initial pressure (Bar)

% Data from Sutton LRE LOX+RP1

Tcc1 = 3571; % Combustion chamber temperature (K)
Tcc2 = 3677; % Combustion chamber temperature (K)

OF1 = 2.24; % OX/FU ratio
OF2 = 2.56; % OX/FU ratio

%% PRELIMINARY NOZZLE
% Choose the maximum length of the additive manufacturing with a certain
% percentage
eps_n = 65; % nozzle area ratio to be decided see document

% Assume isentropic nozzle to found external pressure
fun_Pe = @(Pe) 1/eps_n - ((gamma+1)/2)^(1/(gamma-1)) * (Pe/P_cc_in)^(1/gamma) * sqrt((gamma+1)/(gamma-1) * (1-(Pe/P_cc_in)^((gamma-1)/gamma)));
P_e = fsolve(fun_Pe,0);

% Convert to Pascal
P_e = P_e*10^5;
P_cc_in = P_cc_in*10^5;

% cT and c*
c_t = sqrt(2*gamma^2/(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1))*(1-(P_e/P_cc_in)^((gamma-1)/gamma)))+P_e/P_cc_in*eps_n; %ideal
c_star = sqrt((R/Mm1)*Tcc1/(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1)))); %ideal

% Throat
A_t = T_in/(P_cc_in*c_t);
D_t = sqrt(4*A_t/pi); 

% Exit area
A_e = eps_n*A_t;
D_e = sqrt(4*A_e/pi);


Is = c_t*c_star/g0; % ideal
m_flow = T_in/(Is*g0); % Mass flow rate
v_e = sqrt(2*(gamma/(gamma-1)*R/Mm1*Tcc1*(1-(P_e/P_cc_in)^((gamma-1)/gamma))));
%% CONICAL NOZZLE
alpha = deg2rad(15); % divergent semi-angle
beta = deg2rad(45); % convergent semi-angle

% Mach choosed - I method
% M = 0.2; % desired Mach
% A_cc = A_t/M*((2/(gamma+1))*(1+(gamma-1)/2 * M^2))^((gamma+1)/(2*(gamma-1)));
% D_cc = sqrt(4*A_cc/pi); % combustion chamber diameter

Lstar = 1.143; % characteristic length (m) average 45 inches
Vcc = Lstar * A_t; % combustion chamber volume

% Impose contraction ratio - II method

eps_c = ceil(8*(D_t*100)^(-0.6) + 1.25); % analytic formula from Space propulsion analysis and design
A_cc = eps_c * A_t; % CC area

D_cc = sqrt(4*A_cc/pi);
L_cc = Vcc/A_cc; % CC length


Lc = 1/2*(D_cc-D_t)/tan(beta); % Convergent length
Ld = 1/2*(D_e-D_t)/tan(alpha); % Divergent length

lambda = (1+cos(alpha))/2; % divergence losses coefficient
T_in_2d = T_in*lambda;

%% Viscosity: displacemente thickness (ref: doc displacement thickness on webeep + Non reacting flow slides)
% throatRey         :   Reynolds at throat
% throatCurv        :   curvature at throat
% throatDiameter    :   throat diameter

throatRadius = D_t/2; %throat radius
throatCurv = 1.5*throatRadius;
throatRey = Reynolds(D_t);

modRe = sqrt(throatRadius /throatCurv)*throatRey; % Modified Reynolds

Cd = 1 - ((gamma + 1)/2)^(3/4) * (3.266 - 2.128/(gamma + 1))*modRe^(-1/2) + 0.9428*((gamma - 1)*(gamma + 2)/(gamma + 1)^(1/2))*modRe^(-1);

m_flow_real=m_flow*Cd;
P_a=0;
T_real=m_flow_real*v_e+A_e*(P_e-P_a);
%Cd = mass flow rate real / mass flow rate ideal ~ Area effective / Area throat


%% RAO
Lc_R = Lc;

% MAX PERFORMANCE
Ld_R_maxp = Ld;

% acquire rao angles from graph

teta_fin = deg2rad(4.13);
teta_ini = deg2rad(40.18);

lambda_rao_maxp = 1/2*(1+cos((alpha+teta_fin)/2));
T_in_2d_rao_maxp = T_in*lambda_rao_maxp;

Ln_Rao_maxp = Ld_R_maxp+Lc_R;
Lc = Lc_R;
Ld = Ld_R_maxp;
L_noz = Ln_Rao_maxp;
T_real = T_real*lambda_rao_maxp;



%% MIN LENGTH
% Ld_R_minl = Ld*0.6;
% 
% % acquire rao angles from graph
% 
% teta_fin = deg2rad(12.9);
% teta_ini = deg2rad(38.29);
% 
% alpha_Rao = atan((D_e-D_t)/(2*Ld_R_minl));
% lambda_rao_minl = 1/2*(1+cos((alpha_Rao+teta_fin)/2));
% T_in_2d_rao_minl = T_in*lambda_rao_minl;
% 
% Ln_Rao_minl = Ld_R_minl+Lc_R;

%% COMBUSTION CHAMBER
% Isentropic Nozzle

Is_real = T_real/(m_flow*g0);

L_tot_thrustchamber = L_cc+L_noz; % Engine complete length
m_flow_ox = (OF1/(1+OF1))*m_flow; 
m_flow_f = (1/(1+OF1))*m_flow;


T_tank_ox = 120; % Ox temperature, based on boiling temperature OX at 19 bar (K)
T_tank_fu = 293.15; % Fu temperature based on Standard temperature (K)

rho_ox = 980; % LOX density at 60 bar and T = 120 K (kg/m^3)
rho_fu = 830; % RP-1 density at 60 bar  (kg/m^3)

rho_mix = P_cc_in/(R/Mm1 * Tcc1);

D_cc = sqrt(4*A_cc/pi);
L_cc = Vcc/A_cc;

v_flow = m_flow/(A_cc * rho_mix);
a = sqrt(gamma * R/Mm1 * Tcc1);

M = v_flow/a; % a lot < 0.3 

%% FED SYSTEM
% Distributed pipes pressure drop

L_ox_pipe = 0.83; %[m]
L_fu_pipe = 0.38; %[m] % Consider Regenerative cooling so pipe goes back to front, linear part
roughness_abs = 2e-6; %[m] pipe line roughness of AISI 310 L Stainless Steel
D_pipe_ox = 0.00625; %[m] 1/4 inches
D_pipe_fu = 0.00625; %[m] 1/4 inches
A_pipe_ox = D_pipe_ox^2/4 * pi; %[m^2]
A_pipe_fu = D_pipe_fu^2/4 * pi; %[m^2]
v_ox = m_flow_ox/(rho_ox * A_pipe_ox); %[m/s]
v_fu = m_flow_f/(rho_fu * A_pipe_fu); %[m/s]
mu_ox = 0.0001028; %[Pa * s] % dynamic ox viscosity at tank temperature 120 K
mu_fu = 0.00075; %[Pa * s] % dynamic fu viscosity at tank temperature 293.15 K
Re_ox = rho_ox * v_ox * D_pipe_ox/mu_ox;
Re_fu = rho_fu * v_fu * D_pipe_fu/mu_fu;
roughness_rel_ox = roughness_abs/D_pipe_ox;
roughness_rel_fu= roughness_abs/D_pipe_fu;

% f from Colebrook - White equation

fun_ox = @(f) 1/sqrt(f) + 2 * log10(2.51/(Re_ox * sqrt(f)) + roughness_rel_ox/3.71);
fun_fu = @(f) 1/sqrt(f) + 2 * log10(2.51/(Re_fu * sqrt(f)) + roughness_rel_fu/3.71);

f_fu = fsolve(fun_fu,0.001); % Darcy friction factor
f_ox = fsolve(fun_ox,0.001); % Darcy friction factor

DeltaP_feed_ox = f_ox * L_ox_pipe * rho_ox/2 * v_ox^2/D_pipe_ox;
DeltaP_feed_fu = f_fu * L_fu_pipe * rho_fu/2 * v_fu^2/D_pipe_fu;


% Dynamic pressure drop

DeltaP_dyn_ox = 1/2 * rho_ox * v_ox ^ 2;
DeltaP_dyn_fu = 1/2 * rho_fu * v_fu ^ 2;

% Concentrated Injectors pressure drop
% Choose diameter and number
% Sharp-edged orifice
% triplet impinging stream pattern 2ox 1f
% orifice data
d_inj=0.0014; % (m)
c_d = 0.65;

c_d_ox = c_d;
c_d_fu = c_d;
N_inj_ox = 8;
N_inj_fu = N_inj_ox/2;

k_inj_ox=1/(c_d_ox^2);
k_inj_fu=1/(c_d_fu^2);

A_inj_ox = pi*d_inj^2/4;
A_inj_fu = pi*d_inj^2/4;

% Total area of injection
A_inj_tot_ox = A_inj_ox*N_inj_ox;
A_inj_tot_fu = A_inj_fu*N_inj_fu;


DeltaP_inj_ox = (m_flow_ox/(c_d*A_inj_tot_ox))^2/(2*rho_ox);
DeltaP_inj_fu = (m_flow_f/(c_d*A_inj_tot_fu))^2/(2*rho_fu);

if (DeltaP_inj_fu/P_cc_in)<0.05||(DeltaP_inj_ox/P_cc_in)<0.05
    disp('too low Delta P of injection');
    (DeltaP_inj_ox/P_cc_in)*100;
    (DeltaP_inj_fu/P_cc_in)*100;
end

% Velocity of fuel and oxidier in the injectors
v_ox_inj=sqrt(2*DeltaP_inj_ox/(rho_ox*k_inj_ox));
v_fu_inj=sqrt(2*DeltaP_inj_fu/(rho_fu*k_inj_fu));

%Delta Pressure of the check valve
DeltaP_valve = 2.5 * 10^5; %[Pa]

%Delta Pressure of the entrance in the pipe line
k_entry = 0.5; 
DeltaP_entry_fu = 1/2 * rho_fu * v_fu ^ 2 * k_entry;
DeltaP_entry_ox = 1/2 * rho_ox * v_ox ^ 2 * k_entry;

      
 
%Delta Pressure of the Solenoid valves
K_v = 1.1;
DeltaP_sol_valve_fu = (m_flow_f*15850.3/(K_v * rho_fu))^2 * 6894.76; 
DeltaP_sol_valve_ox = (m_flow_ox*15850.3/(K_v * rho_ox))^2 * 6894.76;

%Delta Pressure loss of the Anular cooling jacket with RP-1
% The hydraulic diameter is chosen...
D_cooling = 0.005; % diameter of cooling channels (m)
D_int = D_e;
D_ext = D_int + D_cooling;
b_a = D_int/D_ext;
% From Table in the reference
xi = 0.667;
D_eff = D_cooling * xi;
% Correct Reynolds and roughness
abs_rough_cool = 50e-6;
roughness_rel_cool= abs_rough_cool/D_eff;
v_cool = m_flow_f/(rho_fu * D_cooling^2 * pi/4);
Re_fu_cool = rho_fu * v_cool * D_eff/mu_fu;

% f from Colebrook - White equation

fun_fu_cool = @(f) 1/sqrt(f) + 2 * log10(2.51/(Re_fu_cool * sqrt(f)) + roughness_rel_cool/3.71);

f_fu_cool = fsolve(fun_fu_cool,0.001); % Darcy friction factor


L_cooling = Ld/cos(alpha)+Lc/cos(beta)+L_cc;
DeltaP_cooling_fu = f_fu_cool * L_cooling * rho_fu/2 * v_fu^2/D_cooling;

P_ox_T = P_cc_in + DeltaP_feed_ox + DeltaP_dyn_ox + DeltaP_inj_ox + DeltaP_valve + DeltaP_entry_ox + DeltaP_sol_valve_ox;
P_fu_T = P_cc_in + DeltaP_feed_fu + DeltaP_dyn_fu + DeltaP_inj_fu + DeltaP_valve + DeltaP_entry_fu + DeltaP_sol_valve_fu+DeltaP_cooling_fu;

%% Tanks
D_cil = 1;
H_cil = 2;

V_cil = D_cil^2/4 * pi * H_cil; % Available volume
V_sys = 0.8 * V_cil; % 80% ov Volume

V_conv = 1/3 * pi * Lc * (D_t^2/4 + D_cc^2/4 + D_t*D_cc/4); % Convergent volume

OFV1=1.59; %volume ox ratio 1
OFV2=1.82; %volume ox ratio 2

P_gasox_fin = P_ox_T - 30e5;
P_gasfu_fin = P_fu_T - 30e5;
V_tot_tank = (V_sys - Vcc - V_conv);
V_tot_tank_ox = (OFV1/(1+OFV1))*V_tot_tank;
V_tot_tank_fu = (1/(1+OFV1))*V_tot_tank;


gamma_gas=1.66; %change when decide gas
V1_gas_fu=V_tot_tank_fu*(P_gasfu_fin/P_fu_T)^(1/gamma_gas);
V_fu=V_tot_tank_fu-V1_gas_fu;

V1_gas_ox=V_tot_tank_ox*(P_gasox_fin/P_ox_T)^(1/gamma_gas);
V_ox=V_tot_tank_ox-V1_gas_ox;

d_tank=0.97*D_cil;
A_tank = pi*d_tank^2/4;
L_tank_ox=V_tot_tank_ox/A_tank;

%% INSULANT
r_i = 0.00625/2;                %[m] Internal tube radius
thick_tube = 2e-3;              % Thickness [m]
r_e = r_i + thick_tube;         %[m] External tube radius

k_t = 15;                       %[W/Km] Tube conductivity
k_is = 0.022;                   %[W/Km] Insulant conductivity 

T_RP1 = 20 + 273.15;            %[K] RP1 temperature
density_RP1 = 830;       %[Kg/m^3] (At T = 298 K)
k_RP1 = 0.115;                   %[W/Km] RP1 conductivity  
cp_RP1 = 1884.1;                %[J/kgK] (At T = 298 K)
dyn_viscosity_RP1 = 0.75e-3;    %[Pa * s] (At T = 298 K)  
L_c_e = 0.97 - r_e*2;  % Characteristic length fuel
                 
t_m = 133;                      %[K] Target's temperature 

T_LOX = 120;                     %[K] LOX temperature
density_LOX = 980;              %[Kg/m^3]
k_LOX = 0.09;                   %[W/Km] LOX conductivity
v_LOX = 7.96;                   %[m/s] LOX velocity
dyn_viscosity_LOX = 0.0001028;           %[Kg/m s]  
cp_LOX = 2302.7;                %[J/kg K] 
Lc_i = r_i * 2;                  %[m]

%LOX REF: https://rocketprops.readthedocs.io/en/latest/lox_prop.html
%RP1 REF: https://kinetics.nist.gov/RealFuels/macccr/macccr2008/Bruno2.pdf
%% Convection coefficient
g = 9.81;
beta = 8.75e-4;
T_is_ext = 17 + 273.15;
kinet_viscosity_RP1 = dyn_viscosity_RP1/density_RP1;
Gr_e = g * beta * (T_RP1-T_is_ext) * L_c_e^3/(kinet_viscosity_RP1^2); % Grashof
Pr_e = (dyn_viscosity_RP1*cp_RP1)/k_RP1;          %Prandtl 
Ra_e = Gr_e * Pr_e;
C_e = 0.023;
a_e = 0.8;
b_e = 0.3;
f_Pr = (1 + (0.49/Pr_e)^(9/16))^(16/9);
Nu_e = 5.75 + 0.75 * (Ra_e/f_Pr)^0.252;   %Nusselt
h_e = (Nu_e*k_RP1)/(L_c_e);                      %External convection coeff --> RP-1 

Rey_i = (density_LOX*Lc_i*v_LOX)/k_LOX;       %Reynolds
Pr_i = (dyn_viscosity_LOX*cp_LOX)/k_LOX;          %Prandtl 
C_i = 0.023;
a_i = 0.8;
b_i = 0.4;
Nu_i = C_i*Rey_i^(a_i)*Pr_i^(b_i);          %Nusselt
h_i = (Nu_i*k_LOX)/Lc_i;                      %Internal convection coeff --> LOX

%% Thickness
K = @(r_is) 2*pi/((1/(h_i*r_i)) + (log(r_e/r_i)/k_t) + (log(r_is/r_e)/k_is) + (1/(h_e*r_is)));
fun = @(r_is) T_LOX + K(r_is)*(T_RP1 - T_LOX)/(2*pi*r_is*h_e) - t_m;
r_is_min = fsolve(fun, 0.005);
thick_is = r_is_min - r_e;
d_pipehole = D_pipe_ox+2*thick_tube+2*thick_is;
A_base_fu = A_tank-d_pipehole^2*pi/4;
L_tank_fu = V_tot_tank_fu/A_base_fu;

L_tot=L_tank_fu+L_tank_ox+L_cc+Lc;

M_f = V_fu*rho_fu/1.02; % Fuel mass
M_ox = V_ox*rho_ox/1.02; % Oxidizer mass

P_control_ox = P_ox_T*2;
P_control_f = P_fu_T*2;

sigma_tank=380e6; %ultimate tensile strenght Al2319 tank (Pa)

thick_tank_ox=P_control_ox*d_tank/(2*sigma_tank);
thick_tank_f=P_control_f*d_tank/(2*sigma_tank);

A_tot_tank_f = A_tank*2 + pi * d_tank * L_tank_fu;
A_tot_tank_ox = A_tank*2 + pi * d_tank * L_tank_ox;

M_fuel_tank = A_tot_tank_f * thick_tank_f * rho_fu;
M_ox_tank = A_tot_tank_ox * thick_tank_ox * rho_ox;

Mm_pg=0.004; % He molar mass (Kg/mol)
T_pg=293.15;

rho_pg_ox=P_ox_T/(R/Mm_pg*T_pg);
rho_pg_fu=P_fu_T/(R/Mm_pg*T_pg);
M_pg_ox=V1_gas_ox*rho_pg_ox;
M_pg_f=V1_gas_fu*rho_pg_fu;

%sigma_tc=41e6;
rho_C_103=8470; % C-103 density (kg/m^3)
%P_control_cc = P_cc_in*2;
%thick_cc=P_control_cc*D_cc/(2*sigma_tc) %thickness thrust chamber
thick_cc=0.003;
V_cc_chamber=D_cc*pi*thick_cc*L_cc;
V_conv = 1/3 * pi * Lc * (D_t^2/4 + D_cc^2/4 + D_t*D_cc/4);
V_div = 1/3 * pi * Ld * (D_t^2/4 + D_e^2/4 + D_t*D_e/4);
V_tot=V_cc_chamber+V_conv+V_div;
M_tc=rho_C_103*V_tot;

%% RAO PLOT


[lengths, nozzle.shapeFun, lossCoefficient, type] = buildNozzle(D_t,D_cc,eps_n);

nozzle.convLenght = lengths(1);
nozzle.divLenght = lengths(2);
nozzle.lenght = lengths(3);
 
% Utility function for better nozzle plotting
nozzle.invShapeFun = @(x) -1 * nozzle.shapeFun(x);

% Plot
fplot(nozzle.shapeFun,[0, nozzle.lenght], 'LineWidth',2, 'Color','k')
hold on;
fplot(nozzle.invShapeFun,[0, nozzle.lenght], 'LineWidth',2, 'Color', 'k')
axis equal
yline(0,'--')

title('Bell Nozzle')



