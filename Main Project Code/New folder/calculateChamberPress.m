function [massFlowRateOx, massFlowRateFu, chamberPress, chamberTemp, chamberGamma, throatGamma, tempFu, tempWall, throatReynolds, hgas] = calculateChamberPress(...
    pressOx, pressFu, estMassFlowRateOx, estMassFlowRateFu, ox, fu, oxTank, fuTank, eFeedLines, eCoolingJacket, injector, engine,...
    tempWall, tempFu, estChamberPress, throatReynolds, chamberGamma, throatGamma)
%calculateChamberPress calculates mass flow rate of fuel and oxidizer and
%chamber pressure for the given tank pressures

densityOx = ox.density;
densityFu = fu.density;
feedLengthOx = oxTank.feedLength;
feedLengthFu = fuTank.feedLength;
feedDiameterOx = oxTank.feedDiameter;
feedDiameterFu = fuTank.feedDiameter;
viscosityOx = ox.viscosity;
viscosityFu = fu.viscosity;
injectorHoleSizeOx = injector.holeSizeOx;
injectorHoleSizeFu = injector.holeSizeFu;
injectorHoleNumberOx = injector.holeNumberOx;
injectorHoleNumberFu = injector.holeNumberFu;
throatArea = engine.throatArea;
initialTankTempOx = oxTank.initialTemp;

injectorAreaOx = pi*injectorHoleSizeOx^2*injectorHoleNumberOx/4;
injectorAreaFu = pi*injectorHoleSizeFu^2*injectorHoleNumberFu/4;
feedAreaOx = pi*feedDiameterOx^2/4;
feedAreaFu = pi*feedDiameterFu^2/4;

Kentrance = 0.5; % From white fluid mechanics for an entry with no bevel
Kvalve = 1.1; % Change this value for a manufacturer provided K

KCdOx = 1/injector.CdOx^2;
KCdFu = 1/injector.CdFu^2;

maxIterations = 150;

i = 0;

edFeedLinesFu = eFeedLines/feedDiameterFu;
edFeedLinesOx = eFeedLines/feedDiameterOx;

while true
    
    vOx = estMassFlowRateOx/(densityOx*feedAreaOx);
    vFu = estMassFlowRateFu/(densityFu*feedAreaFu);

    estReynoldsOx = vOx*densityOx*feedDiameterOx/viscosityOx;
    estReynoldsFu = vFu*densityFu*feedDiameterFu/viscosityFu;
    
    % Colebrook flormula for friction factor
    colebrookOx = @(f) 1/sqrt(f)+2*log10((edFeedLinesOx/3.7)+(2.51)/(estReynoldsOx*sqrt(f)));
    colebrookFu = @(f) 1/sqrt(f)+2*log10((edFeedLinesFu/3.7)+(2.51)/(estReynoldsFu*sqrt(f)));
    
    if estReynoldsOx > 4000 % Turbulent flow
        fOx = fzero(colebrookOx, [0.008, 0.1]);
    elseif estReynoldsOx < 2300 % Laminar flow
        fOx = 64/estReynoldsOx;
    else % Transitional
        fOx = (((estReynoldsOx-2300)/(4000-2300))*(0.1-0.008))+0.008;
    end
    
    if estReynoldsFu > 4000 % Turbulent flow
        fFu = fzero(colebrookFu, [0.008, 0.1]);
    elseif estReynoldsOx < 2300 % Laminar flow
        fFu = 64/estReynoldsFu;
    else % Transitional
        fFu = (((estReynoldsFu-2300)/(4000-2300))*(0.1-0.008))+0.008;
    end
    
    [nozzleFrictionLoss, fNozzle] = calculateNozzleFrictionLosses(engine, fu, estMassFlowRateFu, eCoolingJacket);
    
    checkValveDpOx = calculateCheckValveLoss(estMassFlowRateOx, ox.density);
    checkValveDpFu = calculateCheckValveLoss(estMassFlowRateFu, fu.density);
    
    modRe = sqrt(1/1.5)*throatReynolds;
    
    CdThroat = 1 - ((throatGamma + 1)/2)^(3/4) * (3.266 - 2.128/(throatGamma + 1))*modRe^(-1/2) + 0.9428*((throatGamma - 1)*(throatGamma + 2)/(throatGamma + 1)^(1/2))/modRe;
    
    massFlowRateOx = CdThroat*sqrt(densityOx*(pressOx - estChamberPress - checkValveDpOx)/...
        (KCdOx/(2*injectorAreaOx^2)+(1+Kentrance+Kvalve+fOx*feedLengthOx/feedDiameterOx)/(2*feedAreaOx^2)))/2 + estMassFlowRateOx/2;
    
    massFlowRateFu = CdThroat*sqrt(densityFu*(pressFu - estChamberPress - checkValveDpFu)/...
        (KCdFu/(2*injectorAreaFu^2)+(1+Kentrance+Kvalve+fFu*feedLengthFu/feedDiameterFu)/(2*feedAreaFu^2)+nozzleFrictionLoss))/2 + estMassFlowRateFu/2;
    
    [tempFu, tempWall, hgas] = regenerativeCoolingCalculations(engine, fu, fuTank, oxTank, tempWall, tempFu,massFlowRateOx, massFlowRateFu, estChamberPress, fNozzle);
    
    Output=CEA('problem','rocket','frozen', 'nfr', 1,'o/f',massFlowRateOx/massFlowRateFu,'sup,ae/at',engine.expansionRatio, 'p(bar)',estChamberPress/10^5,...
        'reactants','fuel','RP-1(L)','C',1.,'H',1.9423,'t(k)',tempFu,'oxid','O2(L)','O',2,...
        't(k)',initialTankTempOx,'output','default','transport','end');
    
    chamberGamma = Output.output.froz.gamma(2);
    throatGamma = Output.output.froz.gamma(2);
    chamberTemp = Output.output.froz.temperature(1);
    R = 8.314*1000/Output.output.froz.mw(1);
    throatReynolds = Output.output.froz.cstar(2)*engine.throatDiameter*Output.output.froz.density(2)/(Output.output.froz.viscosity(2)*0.0001);
    
    chamberPress = (massFlowRateFu + massFlowRateOx)*sqrt(R*chamberTemp)/...
        (throatArea*sqrt(chamberGamma))*(((chamberGamma+1)/2)^((chamberGamma+1)/(2*(chamberGamma-1))))/2 + estChamberPress/2;
    
    % Check convergence
    if  abs((massFlowRateOx - estMassFlowRateOx)/massFlowRateOx) < 0.001 && ...
        abs((massFlowRateFu - estMassFlowRateFu)/massFlowRateFu) < 0.001 && ...
        abs((chamberPress - estChamberPress)/chamberPress) < 0.001
        % Add warning if reynolds belongs to [2000; 4000]
        %fprintf('Iterations required: %d.\n',i);
        break;
    end
    
    i = i + 1;
    if i == maxIterations
        fprintf('Warning: Max iterations reached.\n');
        break
    end
    
    estMassFlowRateOx = massFlowRateOx;
    estMassFlowRateFu = massFlowRateFu;
    estChamberPress = chamberPress;
end