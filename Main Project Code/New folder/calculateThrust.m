function thrust = calculateThrust(chamberPress, throatArea, gamma, expansionRatio, teta_fin, alpha)
%calculateThrust calculates the produced thrust

fun = @(Pe) 1/expansionRatio - ((gamma+1)/2)^(1/(gamma-1)) * (Pe/chamberPress)^(1/gamma) * sqrt((gamma+1)/(gamma-1) * (1-(Pe/chamberPress)^((gamma-1)/gamma)));

exitPress = fzero(fun, [0, chamberPress/expansionRatio]);

Ct = sqrt(2*gamma^2/(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1))*(1-(exitPress/chamberPress)^((gamma-1)/gamma)))+exitPress/chamberPress*expansionRatio;

%% RAO (ref: exercise session 03/04/2024)
% acquire rao angle from graph --> WE HAVE TO CHANGE IT (ask to Pit)

lambdaRAO = 1/2*(1+cosd((alpha+teta_fin)/2)); %Losses
%% Thrust

thrust = chamberPress*throatArea*Ct*lambdaRAO;

end