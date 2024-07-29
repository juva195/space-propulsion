function [lengths, shapeFun, lossCoefficient, type] = buildNozzle(dThroat,dConv,exRatio,type)
% FUNCTION NAME: 
% buildNozzle
%
% DESCRIPTION:
%   The function computes the nozzle shape and 2D losses for the selected
%   nozzle kind
%
% INPUT:
%   1) dThroat - (1x1) Throat diameter in [m]
%   2) dConv - (1x1) Convergent section diameter in [m]
%   3) exRatio - (1x1) Expansion ratio, defined as areaExit/areaThroat
%   4) kind - (1x1 string) Optional nozzle kind specification, if not
%   present the function computes the "conic" nozzle
%
% OUTPUT:
%   1) lengths - (1x3) Convergent, divergent and total (convergent + divergent sections) lenghts in [m]
%   2) shapeFun - (1x1 function) Piecewise function that computes the y
%   coordinate of the nozzle external wall when given axial coordinate 
%   3) lossCoefficient - (1x1) 2D loss coefficient 
%
% EXAMPLE:
%   [lenght, shapeFun, lossCoefficient] = buildNozzle(dThroat,dChamber,exRatio,"conic");
%   [~, shapeFun] = buildNozzle(dThroat,dChamber,exRatio)
%
%   To plot the nozzle shape:
%   fplot(shapeFun, [0, lenght])
%   yline(0,'--')
%   axis equal
%   axis padded
%
% CALLED FUNCTIONS:
%   raoNozzle
%
% REVISION HISTORY:
%   15/12/2023 - Francesco Miccoli
%       * First implementation
%
% LICENSE:
%   Copyright © 2023, Skyward Experimental Rocketry, PRP department
%   All rights reserved
%   SPDX-License-Identifier: GPL-3.0-or-later
%
% CONTACT:
%   Skyward Experimental Rocketry: info@skywarder.eu
%   Francesco Miccoli: francesco.miccoli@skywarder.eu
%
%
% ASSUMPTIONS AND LIMITATIONS:
%   None
%
% REFERENCES:
%   None
%

%% INPUT CHECK

narginchk(3,4)
if nargin == 3
    type = "rao";
end

%% BUILD NOZZLE

% Impose conic nozzle for certain expansion ratios
if exRatio < 4.5 && type ~= "conic"
    type = "conic";
    warning("The nozzle type has been automatically set to 'conic' because of it's low expansion ratio")
end

switch type

    case "conic"
        
        % Set convergent and divergent sections angles
        convAngle = deg2rad(45);
        divAngle = deg2rad(20);
        
        % Compute exit station diameter
        dExit = sqrt( exRatio * dThroat^2 );
        
        % Compute lenghts
        convLenght = (dConv-dThroat)/2 /sin(convAngle) * cos(convAngle);
        divLenght = (dExit-dThroat)/2 /sin(divAngle) * cos(divAngle);
        lenght = convLenght + divLenght;
        lengths = [convLenght, divLenght, lenght];
        
        % Build piecewise function
        shapeFun = @(x) (dConv/2 - x./cos(convAngle) .* sin(convAngle)) .* ((0<=x) & (x<=convLenght)) + ...
                        (dThroat/2 + (x-convLenght)./cos(divAngle) .* sin(divAngle)) .* ((convLenght<x) & (x<=lenght));

        % Compute 2D losses
        lossCoefficient = (1 + cos(divAngle)) /2;

    case "rao"

        [lengths,shapeFun,thetaE] = raoNozzle(dThroat,dConv,exRatio,80);

        % Compute 2D losses
        lossCoefficient = (1 + cos(thetaE)) /2;

    otherwise

        error('Please either omit "type" input or use a supported "type", for more see function help')

end


end