function checkValveDp = calculateCheckValveLoss(massFlowRate, density)
%calculateCheckValveLoss calculates de pressure drop in the check valve

flowRates = [0 5 6.217949 7.136752 7.927350 8.482906 9.914530 11.538462 12.991453 14.957265]; % [L/min]

pressureDrops = [0.5 0.5 0.545455 0.619318 0.744318 0.886364 1.318182 1.789773 2.221591 2.784091]*10^5; % Pa

flowRate = massFlowRate/density*1000*60;

checkValveDp = interp1(flowRates, pressureDrops, flowRate, 'linear','extrap');

if isnan(checkValveDp)
    error('Flow rate higher than table.');
end

end