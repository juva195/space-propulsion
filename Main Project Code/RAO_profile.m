function [xPoints, yPoints, totalLength, length, alpha] = RAO_profile(engine)
%RAO_profile calculates the profile of a RAO nozzle

Rt = engine.throatDiameter/2;
expansionRatio = engine.expansionRatio;
sigma_n = engine.sigma_n;
sigma_e = engine.sigma_e;
nozzleLength = engine.nozzleLength;
nPoints = engine.nPoints;
combustionChamberLength = engine.combustionChamberLength;
combustionChamberDiameter = engine.combustionChamberDiameter;

Ln = nozzleLength * Rt * (sqrt(expansionRatio) - 1) / tand(15);

initialAngle = -45;

xMin = 1.5*Rt*cosd(initialAngle);
yMin = 1.5*Rt*sind(initialAngle)+1.5*Rt+Rt;

length45 = combustionChamberDiameter/2-yMin;

totalLength = Ln + xMin + combustionChamberLength + length45;

Nx = 0.382*Rt*cosd(sigma_n-90);
Ny = 0.382*Rt*sind(sigma_n-90)+0.382*Rt+Rt;

Ex = Ln;
Ey = sqrt(expansionRatio)*Rt;

m1 = tand(sigma_n);
m2 = tand(sigma_e);
C1 = Ny - m1*Nx;
C2 = Ey - m2*Ex;

Qx = (C2-C1)/(m1-m2);
Qy = (m1*C2-m2*C1)/(m1-m2);

nPointsC = round(nPoints * 0.15);

nPointsN = round(nPoints * 0.4);

nPointsT = nPoints - nPointsC - nPointsN;

coolingChannelNumber = [ones([1, nPointsC]) * engine.coolingChannelNumberC, ones([1,nPointsT]) * engine.coolingChannelNumberT, ones([1,nPointsN-1]) * engine.coolingChannelNumberN];

channelThickness = [ones([1, nPointsC]) * engine.channelThicknessC, ones([1,nPointsT]) * engine.channelThicknessT, ones([1,nPointsN-1]) * engine.channelThicknessN];

xPointsC = linspace(-xMin - combustionChamberLength, -xMin - length45, nPointsC);

xPointsT = linspace(-xMin - length45, Nx, nPointsT + 2);

xPointsN = linspace(Nx, Ln, nPointsN);

xPoints = [xPointsC, xPointsT(2:end-1), xPointsN];
yPoints = [];

length = xPoints(2:end) - xPoints(1:end-1);

for x = xPoints
    if x <= -xMin-length45
        yPoints(end+1) = combustionChamberDiameter/2;
    elseif x< -xMin
        yPoints(end+1) = combustionChamberDiameter/2 - (x + xMin + length45);
    elseif x <= 0
        sigma = -acosd(x/(1.5*Rt));
        yPoints(end+1) = 1.5*Rt*sind(sigma)+1.5*Rt+Rt;
    elseif x <= Nx
        sigma = -acosd(x/(0.382*Rt));
        yPoints(end+1) = 0.382*Rt*sind(sigma)+0.382*Rt+Rt;
    else
        fun = @(t) (1-t^2)*Nx+2*(1-t)*t*Qx+t^2*Ex - x;
        tPoint = fzero(fun, [0,1]);
        yPoints(end+1) = (1-tPoint^2)*Ny+2*(1-tPoint)*tPoint*Qy+tPoint^2*Ey;
    end
end

xPoints = xPoints + xMin + combustionChamberLength + length45;

alpha = atand((Ey-Rt)/Ln);

end