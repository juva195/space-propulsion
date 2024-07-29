clear 
close all
clc

%% INPUT DATA
M_ox=2.319737663569266e+02;
M_f=1.487838923032807e+02;
M_tc_100=6.804413421818744; 
M_p=M_ox+M_f;
MRpt= 0.85; % mass ratio propellant over total for upper stages
Mtot=M_p/MRpt;
Mfin=Mtot-M_p;

R = 8.314462618; % universal gas costant (J/(K*mol))
g0 = 9.81; % gravity acceleration (m/s^2)
T_in = 1000; % Initial Thrust (N)

% Data from Sutton LRE LOX+RP1

Mm1 = 21.9e-3; % Molar mass (Kg/mol)
Mm2 = 23.3e-3; % Molar mass (Kg/mol)

gamma = 1.24; % Specific heat ratio of the mixture

Tcc1 = 3571; % Combustion chamber temperature (K)
Tcc2 = 3677; % Combustion chamber temperature (K)

P_cc_in = 50; %Combustion chamber initial pressure (Bar)

OF1 = 2.24; % OX/FU ratio
OF2 = 2.56; % OX/FU ratio

%% PRELIMINARY NOZZLE

Is=zeros(451,1);
L_tc=zeros(451,1);
DM=zeros(451,1);
DV=zeros(451,1);

for eps_n=50:500
    % assumption isentropic nozzle
    R = 8.314462618; % universal gas costant (J/(K*mol))
    g0 = 9.81; % gravity acceleration (m/s^2)
    T_in = 1000; % Initial Thrust (N)
    
    % Data from Sutton LRE LOX+RP1
    
    Mm1 = 21.9e-3; % Molar mass (Kg/mol)
    Mm2 = 23.3e-3; % Molar mass (Kg/mol)
    
    gamma = 1.24; % Specific heat ratio of the mixture
    
    Tcc1 = 3571; % Combustion chamber temperature (K)
    Tcc2 = 3677; % Combustion chamber temperature (K)
    
    P_cc_in = 50; %Combustion chamber initial pressure (Bar)
    
    OF1 = 2.24; % OX/FU ratio
    OF2 = 2.56; % OX/FU ratio

    fun_Pe = @(Pe) 1/eps_n - ((gamma+1)/2)^(1/(gamma-1)) * (Pe/P_cc_in)^(1/gamma) * sqrt((gamma+1)/(gamma-1) * (1-(Pe/P_cc_in)^((gamma-1)/gamma)));
    
    P_e = fsolve(fun_Pe,0);
    
    % conversion to Pascal
    P_e = P_e*10^5;
    P_cc_in = P_cc_in*10^5;
    
    c_t = sqrt(2*gamma^2/(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1))*(1-(P_e/P_cc_in)^((gamma-1)/gamma)))+P_e/P_cc_in*eps_n; %ideal
    c_star = sqrt((R/Mm1)*Tcc1/(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1)))); %ideal
    
    A_t = T_in/(P_cc_in*c_t); % Throat area
    D_t = sqrt(4*A_t/pi);
    
    A_e = eps_n*A_t; % Exit area
    D_e = sqrt(4*A_e/pi);
   
    alpha = deg2rad(15); % divergent semi angle
    beta = deg2rad(45); % convergent semi angle
    
    % M = 0.2; % desired Mach
    % A_cc = A_t/M*((2/(gamma+1))*(1+(gamma-1)/2 * M^2))^((gamma+1)/(2*(gamma-1)));
    % D_cc = sqrt(4*A_cc/pi); % combustion chamber diameter
    
    Lstar = 1.143; % characteristic length (m) (40/50 inches)
    Vcc = Lstar*A_t; % combustion chamber volume
    
    % Impose contraction ratio - New method
    eps_c = 8*D_t^(-0.6)+1.25; % analytic formula from Space propulsion analysis and design
    A_cc = eps_c * A_t; % CC area
    
    D_cc = sqrt(4*A_cc/pi);
    L_cc = Vcc/A_cc; % CC length
    
    
    Lc = 1/2*(D_cc-D_t)/tan(beta); % Convergent length
    Ld = 1/2*(D_e-D_t)/tan(alpha); % Divergent length
    L_noz = Lc+Ld;
    L_tc(eps_n-49)=L_noz+L_cc;
    
    lambda = (1+cos(alpha))/2; % divergence losses coefficient
    T_in_2d = T_in*lambda;
      
    
    v_e = sqrt(2*(gamma/(gamma-1)*R/Mm1*Tcc1*(1-(P_e/P_cc_in)^((gamma-1)/gamma))));
    
    Is(eps_n-49) = c_t*c_star/g0;

    rho_C_103=8470; % C-103 density (kg/m^3)
    thick=0.003; %thickness thrust chamber
    V_cc=D_cc*pi*thick*L_cc;
    V_conv = 1/3 * pi * Lc * (D_t^2/4 + D_cc^2/4 + D_t*D_cc/4);
    V_div = 1/3 * pi * Ld * (D_t^2/4 + D_e^2/4 + D_t*D_e/4);
    V_tot=V_cc+V_conv+V_div;
    M_tc=rho_C_103*V_tot;
    DM(eps_n-49)=M_tc-M_tc_100;
    DV(eps_n-49)=Is(eps_n-49)*g0*log((DM(eps_n-49)+Mfin+M_p)/(Mfin+DM(eps_n-49)));
end
eps_n=50:500;
plot(eps_n,DV)
figure()
plot(L_tc,eps_n)
perc_DV=(DV-DV(50))./DV(50);
perc_Is=(Is-Is(50))./Is(50);
plot(eps_n,perc_DV)
hold on
perc_M=DM./Mfin;
plot(eps_n,perc_M)
legend('perc DV','perc M')
figure()
plot(eps_n,perc_Is-perc_M)
legend('perc DV-perc M')
[M,I]=max(perc_DV);
[M,I]=max(perc_Is)
eps_n(I)
plot(eps_n,L_tc)



