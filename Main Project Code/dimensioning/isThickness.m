%% Isolant's thickness 
clear
clc
close all

%% Data 
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
density_LOX = 850;              %[Kg/m^3]
k_LOX = 0.09;                   %[W/Km] LOX conductivity
v_LOX = 7.96;                   %[m/s] LOX velocity
dyn_viscosity_LOX = 7e-5;           %[Kg/m s]  
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

%% Thickness Yesterday
K = @(r_is) 2*pi/((1/(h_i*r_i)) + (log(r_e/r_i)/k_t) + (log(r_is/r_e)/k_is) + (1/(h_e*r_is)));
fun = @(r_is) T_LOX + K(r_is)*(T_RP1 - T_LOX)/(2*pi*r_is*h_e) - t_m;
r_is_min = fsolve(fun, 0.005);
thick_is = r_is_min - r_e;

%% Thickness Today
% L_tank = 0.6758;
% S_i = 2*pi*r_i*L_tank;
% S_e = 2*pi*r_e*L_tank;
% 
% R_kt = 1/(2*pi*k_t*L_tank)*log(r_e/r_i);
% R_k_is = @(thick_is) 1/(2*pi*k_is*L_tank)*log((r_e + thick_is)/r_e);
% R_convLOX = 1/(h_i*S_i);
% R_convRP1 = 1/(h_e*S_e);
% 
% T_LOX_initial = 133;
% mass_OX = 231.27;
% Q_LOX = mass_OX*cp_LOX*(T_LOX_initial - T_LOX);
% %%
% Q = @(thick_is) (T_RP1 - T_LOX)/(R_convRP1 + R_convLOX + R_k_is(thick_is) + R_kt) - Q_LOX;
% thick_is_min = fsolve(Q, 0.001);
% % 
% 
% 
% 
% 
