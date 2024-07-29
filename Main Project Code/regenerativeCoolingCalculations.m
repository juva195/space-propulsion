function [tempFu, tempWall, hgas] = regenerativeCoolingCalculations(engine, fu, fuTank, oxTank, estTempWall,...
    estTempFu,massFlowRateOx, massFlowRateFu, chamberPress, fNozzle)
%regenerativeCoolingCalculations calculates the temperature of the fuel
%after the regenerative cooling

yPoints = engine.yPoints;
xPoints = engine.xPoints;
throatArea = engine.throatArea;
Dstar = engine.throatDiameter;
nPoints = engine.nPoints;
wallThickness = engine.coolingWallThickness;
channelThickness = engine.channelThickness;
conductivity = engine.materialConductivity;
coolingDistance = engine.coolingDistance;
channelNumberT = engine.coolingChannelNumberT;
channelNumberN = engine.coolingChannelNumberN;
channelNumberC = engine.coolingChannelNumberC;
channelThicknessT = engine.channelThicknessT;
channelThicknessN = engine.channelThicknessN;
channelThicknessC = engine.channelThicknessC;

throatBeginIndex = find(xPoints < engine.throatBegin,1, 'last');
throatEndIndex = find(xPoints > engine.throatEnd,1, 'first');
channelNumber(1:throatBeginIndex-1) = channelNumberC;
channelNumber(throatBeginIndex:throatEndIndex-1) = channelNumberT;
channelNumber(throatEndIndex:nPoints-1) = channelNumberN;
channelThickness(1:throatBeginIndex-1) = channelThicknessC;
channelThickness(throatBeginIndex:throatEndIndex-1) = channelThicknessT;
channelThickness(throatEndIndex:nPoints-1) = channelThicknessN;

diameter = [yPoints(1:end-1)*2];
area = pi*diameter.^2/4;
length = engine.length;
surfaceAreaInt = pi*diameter.*length./channelNumber;
surfaceAreaExt = pi*(diameter+2*wallThickness).*length./channelNumber;
Dh = 2*channelThickness;

maxIterations = 10;
i = 0;

firstCoolingIndex = find(xPoints(xPoints <= coolingDistance),1,'last');

while true
    
    tempFu = fuTank.initialTemp;

    Output=CEA('problem','rocket','frozen', 'nfr', 1,'o/f',massFlowRateOx/massFlowRateFu,'sup,ae/at',engine.expansionRatio, 'p(bar)',chamberPress/10^5,...
        'reactants','fuel','RP-1(L)','C',1.,'H',1.9423,'t(k)',estTempFu,'oxid','O2(L)','O',2,...
        't(k)',oxTank.initialTemp,'output','default','transport','end');
    
    gamma = Output.output.froz.gamma(2);
    chamberTemp = Output.output.froz.temperature(1);
    cstar = Output.output.froz.cstar(2);
    
    machFun = @(m) 1./m .* (2/(gamma+1) .* (1 + (gamma-1)/2 .* m.^2)) .^ ((gamma+1)/(2*(gamma-1))) - area/throatArea;
    
    options = optimoptions('fsolve','Display','off');
    
    guess = (xPoints >= 0)*2 + (xPoints < 0)*0.5;

    mach = fsolve(machFun, guess(1:end-1), options);
    
    tempGas = (1 + (gamma-1)/2 * mach.^2).^-1 * chamberTemp;
    
    pressGas = (1 + (gamma-1)/2 * mach.^2).^(-gamma/(gamma-1)) * chamberPress;
    
    tempStar = tempGas .* (1+0.032*mach.^2+0.58*(estTempWall./tempGas-1));
    
    Pr = zeros([1,nPoints-1]);
    Cp = zeros([1,nPoints-1]);
    Mu = zeros([1,nPoints-1]);
    
    parfor k = 1:nPoints-1
        Output = 0;
        Output=CEA('problem','tp','o/f',massFlowRateOx/massFlowRateFu,'case','fastrak(O2/RP-1)', 'p(bar)',pressGas(k)/10^5,...
            'T(K)',tempStar(k),'reactants','fuel','RP-1(L)','C',1.,'H',1.9423,'wt%',100., 'h,cal/mol',-5430.,'t(k)',estTempFu,'oxid','O2(L)','O',2,...
            'wt%',100., 'h,cal/mol',-3032.,'t(k)',oxTank.initialTemp,'output','default','transport','end');
        Pr(k) = Output.output.prandtl.eql;
        Cp(k) = Output.output.cp*1000;
        Mu(k) = Output.output.viscosity*0.0001;
    end
    
    sigma = ((1/2.*estTempWall./chamberTemp.*(1+(gamma-1)/2.*mach.^2)+1/2).^0.68.*(1+(gamma-1)/2.*mach.^2).^0.12).^-1;
    
    hgas = 0.026/Dstar^0.2*(Mu.^0.2.*Cp./Pr.^0.6).*(chamberPress/cstar).^0.8.*(1/1.5).^0.1.*(throatArea./area).^0.9.*sigma;
    
    Vfu = channelNumber.*4.*massFlowRateFu./(fu.density*pi*((diameter+wallThickness*2+channelThickness*2).^2-(diameter+wallThickness*2).^2));
    
    Re = Vfu.*fu.density.*Dh/fu.viscosity;
    
    DiDo = (diameter+wallThickness*2)./(diameter+wallThickness*2+channelThickness*2);
    
    NuiTable = [11.56 7.37 5.74 4.86];
    
    DiDoTable = [0.1 0.25 0.5 1];
    
    Nui = interp1(DiDoTable, NuiTable, DiDo);
    
    K1 = 1 + 3.4*fNozzle;
    
    K2 = 11.7 + 1.8/fu.prandtl^(1/3);
    
    Nut = (fNozzle/8).*Re*fu.prandtl./(K1+K2*(fNozzle/8).^0.5*(fu.prandtl^(2/3)-1));
    
    Nu = zeros([1, nPoints-1]);
    
    for j = 1:(nPoints-1)
        if Re(j) < 2300
            Nu(j) = Nui(j);
        elseif Re(j) > 4000
            Nu(j) = Nut(j);
        else
            NuMin = Nui(j);
            K1j = 1 + 3.4*fNozzle(j);
            NuMax = (fNozzle(j)/8).*4000*fu.prandtl./(K1j+K2*(fNozzle(j)/8).^0.5*(fu.prandtl^(2/3)-1));
            Nu(j) = interp1([2300 4000], [NuMin NuMax], Re(j));
        end
    end
    hfu = fu.k./(2*channelThickness).*(Nu);
    
    UA = 1./(1./(hgas.*surfaceAreaInt)+1./(hfu.*surfaceAreaExt)+wallThickness./((surfaceAreaExt+surfaceAreaInt)/2*conductivity));
    
    q = zeros([1,nPoints-1]);

    for l = 1:channelNumberT
        if mod(l,2)
            order = throatEndIndex:-1:throatBeginIndex+1;
        else
            order = throatBeginIndex+1:throatEndIndex;
        end
        for k = order
            LMTD = @(To) (To - tempFu)/(log(tempGas(k-1) - tempFu) - log(tempGas(k-1) - To));
            eq = @(To) UA(k-1)*LMTD(To) - massFlowRateFu*fu.Cp*(To-tempFu);
            nextTempFu = fzero(eq, [tempFu+0.00001, tempGas(k-1)]);
            q(k-1) = q(k-1) + massFlowRateFu*fu.Cp*(nextTempFu-tempFu);
            tempFu = nextTempFu;
        end
    end

    for l = 1:channelNumberN
        if mod(l,2)
            order = firstCoolingIndex:-1:throatEndIndex+1;
        else
            order = throatEndIndex+1:firstCoolingIndex;
        end
        for k = order
            LMTD = @(To) (To - tempFu)/(log(tempGas(k-1) - tempFu) - log(tempGas(k-1) - To));
            eq = @(To) UA(k-1)*LMTD(To) - massFlowRateFu*fu.Cp*(To-tempFu);
            nextTempFu = fzero(eq, [tempFu-1, tempGas(k-1)]);
            q(k-1) = q(k-1) + massFlowRateFu*fu.Cp*(nextTempFu-tempFu);
            tempFu = nextTempFu;
        end
    end

    for l = 1:channelNumberC
        if mod(l,2)
            order = throatBeginIndex:-1:2;
        else
            order = 2:throatBeginIndex;
        end
        for k = order
            LMTD = @(To) (To - tempFu)/(log(tempGas(k-1) - tempFu) - log(tempGas(k-1) - To));
            eq = @(To) UA(k-1)*LMTD(To) - massFlowRateFu*fu.Cp*(To-tempFu);
            nextTempFu = fzero(eq, [tempFu-1, tempGas(k-1)]);
            q(k-1) = q(k-1) + massFlowRateFu*fu.Cp*(nextTempFu-tempFu);
            tempFu = nextTempFu;
        end
    end
    
    tempWall = tempGas - q./(surfaceAreaInt.*channelNumber.*hgas);

    if abs((tempFu-estTempFu)/tempFu) < 0.01
        break;
    end
    
    estTempWall = tempWall;
    estTempFu = tempFu;
    
    %Output=CEA('problem','rocket','equilibrium','o/f',massFlowRateOx/massFlowRateFu,'case','fastrak(O2/RP-1)', 'p(bar)',chamberPress/10^5,...
        %'sub',[area(xPoints <= 0)/throatArea],'sup',[area(xPoints > 0)/throatArea],'reactants','fuel','RP-1(L)','C',1.,'H',1.9423,'wt%',100., 'h,cal/mol',-5430.,'t(k)',estTempFu,...
        %'oxid','O2(L)','O',2,'wt%',100., 'h,cal/mol',-3032.,'t(k)',initialTankTempOx,'output','calories','transport','end');
    
    i = i + 1;
    if i >= maxIterations
        fprintf('Warning: Max cooling iterations reached.');
        break;
    end
end
end