function [nozzleFrictionLoss, fNozzle] = calculateNozzleFrictionLosses(engine, fu, massFlowRateFu, eCoolingJacket)
%calculateNozzleFrictionLosses calculates the friction losses in the nozzle
%cooling jacket

xPoints = engine.xPoints;
yPoints = engine.yPoints;
totalLength = engine.totalLength;
nPoints = engine.nPoints;
wallThickness = engine.coolingWallThickness;
diameter = [yPoints(1:end-1)*2];
length = engine.length;

coolingDistance = engine.coolingDistance;
channelNumberT = engine.coolingChannelNumberT;
channelNumberN = engine.coolingChannelNumberN;
channelNumberC = engine.coolingChannelNumberC;
channelThicknessT = engine.channelThicknessT;
channelThicknessN = engine.channelThicknessN;
channelThicknessC = engine.channelThicknessC;


throatBeginIndex = find(xPoints < engine.throatBegin,1, 'last');
throatEndIndex = find(xPoints > engine.throatEnd,1, 'first');
channelNumber(1:throatBeginIndex) = channelNumberC;
channelNumber(throatBeginIndex+1:throatEndIndex-1) = channelNumberT;
channelNumber(throatEndIndex:nPoints-1) = channelNumberN;
channelThickness(1:throatBeginIndex) = channelThicknessC;
channelThickness(throatBeginIndex+1:throatEndIndex-1) = channelThicknessT;
channelThickness(throatEndIndex:nPoints-1) = channelThicknessN;

firstCoolingIndex = find(xPoints(xPoints <= coolingDistance),1,'last');

area = pi*((diameter+wallThickness*2+channelThickness*2).^2-(diameter+wallThickness*2).^2)/4;

Vfu = channelNumber.*massFlowRateFu./(fu.density.*area);

diamRatio = (diameter + 2*wallThickness)./(diameter + 2*wallThickness + 2*channelThickness);

correctionFactorTable = [0.742, 0.716, 0.693, 0.676, 0.670, 0.667, 0.667];
diamRatios = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1];

correctionFactor = interp1(diamRatios, correctionFactorTable, diamRatio);

Dh = 2*channelThickness;

Deff = Dh.*correctionFactor;

ReynoldsNumbers = Vfu .* Deff * fu.density / fu.viscosity;

edCoolingJacket = eCoolingJacket./Deff;

fNozzle = zeros([1, nPoints-1]);

%options = optimoptions('fsolve','Display','off');

for i = 1:firstCoolingIndex
    Re = ReynoldsNumbers(i);
    colebrook = @(f) 1./sqrt(f)+2*log10((edCoolingJacket(i)./3.7)+(2.51)./(Re.*sqrt(f)));
    if Re < 2300
        fNozzle(i) = 64/Re;
    elseif Re > 4000
        fNozzle(i) = fzero(colebrook, [0.008, 0.1]);
    else
        colebrook = @(f) 1./sqrt(f)+2*log10((edCoolingJacket(i)./3.7)+(2.51)./(4000.*sqrt(f)));
        fMax = fzero(colebrook, [0.008, 0.1]);
        Rey = [2300 4000];
        f = [64/Re fMax];
        fNozzle(i) = interp1(Rey, f, Re);
    end
end
    
nozzleFrictionLoss = sum(fNozzle.*length./(2*Dh.*area.^2));

end