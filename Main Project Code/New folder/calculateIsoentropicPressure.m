function press = calculateIsoentropicPressure(initialMass, initialTankPress, mass, density, gamma, Vt)
%calculateIsoentropicPressure calculates the pressure in a tank with
%isoentropic expansion

Vi = initialMass/density;
V = mass/density;

press = initialTankPress*((Vt-Vi)/(Vt-V))^gamma;

end