figure();
plot(time(1:end-1), thrust(2:end)); % Thrust

figure();
plot(time(1:end-1), chamberPress(2:end)); % Chamber Pressure

figure();
plot(time(1:end-1), massFlowRateOx(2:end)./massFlowRateFu(2:end)); % O/F ratio
yline(2.24);

figure(); % Pressure drops
plot(time(1:end-1), (pressOx(2:end)-chamberPress(2:end))/10^5, 'b');
hold on
plot(time(1:end-1), (pressFu(2:end)-chamberPress(2:end))/10^5, 'r');

figure();
plot(time(1:end-1), thrust(2:end)./(9.81*(massFlowRateOx(2:end)+massFlowRateFu(2:end)))); % ISP

figure();
plot(time(1:end-1) ,tempFu(2:end)); % Fuel temperature

figure();
plot(engine.xPoints(1:end-1), tempWall(end,:)); % Wall temperature on last time step
xline(engine.throatBegin);
xline(engine.coolingDistance);
xline(engine.throatEnd);

figure();
plot(time(1:end-1), max(tempWall(2:end,:),[],2)); % Max wall temp on time

figure();
plot(time(1:end-1), chamberTemp(2:end)); % Chamber temperature