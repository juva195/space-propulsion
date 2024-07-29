function [lengths,shapeFun,exitAngle] = raoNozzle(dThroat,dConv,exRatio,relativeLenght)
% FUNCTION NAME: 
% raoNozzle
%
% DESCRIPTION:
%   The function computes the Rao bell nozzle approximation for the given
%   data. The function should not be called directly, but by buildNozzle.m
%
% INPUT:
%   1) dThroat - (1x1) Throat diameter in [m]
%   2) dConv - (1x1) Convergent section diameter in [m]
%   3) exRatio - (1x1) Expansion ratio, defined as areaExit/areaThroat
%   4) relativeLenght - (1x1) Value representing the percentage of lenght
%   with respect to a 15deg conic length, either 60,70,80,90,100
%
% OUTPUT:
%   1) lengths - (1x3) Convergent, divergent and total (convergent + divergent sections) lenghts in [m]
%   2) shapeFun - (1x1 function) Piecewise function that computes the y
%   coordinate of the nozzle external wall when given axial coordinate 
%   3) exitAngle - (1x1) Exit angle of parabola in [rad]
%
% EXAMPLE:
%   [lenght, shapeFun, thetaE] = buildNozzle(dThroat,dChamber,exRatio,80);
%
% CALLED FUNCTIONS:
%   None
%
% REVISION HISTORY:
%   26/02/2024 - Francesco Miccoli
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

narginchk(4,5);
if ~ismember(relativeLenght,[70,80,90,100])
    error("Please select a valid relative lenght!")
end

%% PROBLEM DATA

dExit = sqrt(exRatio) * dThroat;

%% CONVERGENT SECTION - THROAT ENTRANCE 

% Shape fun - a circle arc with radius 1.5 * R_t
yConvRoundFun = @(x) dThroat/2 * (1.5 + 1) - sqrt((dThroat * 3/4)^2 - (x).^2);

%% CONVERGENT SECTION - STRAIGHT PART 

% Converngent angle
convAngle = deg2rad(45);

% Lenght (on the x axis) of the round part of the convergent section
convRoundLenght = dThroat * 3/4 * sin(convAngle);

% Y coordinate at the straight-round interface of the convergent section
yConvRound = yConvRoundFun(convRoundLenght);

% Straight part shape fun and lenght
convStraightFun = @(x) yConvRound + x./cos(convAngle) .* sin(convAngle);
convStraightLength = (dConv/2-yConvRound)/asin(convAngle) * cos(convAngle);

% Complete convergent part shape fun and lenght
convLength = convRoundLenght + convStraightLength;
convFun = @(x) convStraightFun(convLength-x-convRoundLenght) * ((x>=0) & (x<=convStraightLength)) + yConvRoundFun(convLength-x) * ((x>convStraightLength) & (x<=convLength));

%% DIVERGENT SECTION - THROAT EXIT

% Shape fun - a circle arc with radius 0.382 * R_t
divRoundFun = @(x) dThroat/2 * (0.382 + 1) - sqrt((dThroat/2 * 0.382)^2 - (x-convLength).^2);

%% DIVERGENT SECTION - BELL

% Load Rao Angles interpolated data
load("thetaE.mat","thetaE")
load("thetaN.mat","thetaN")

% Check if thetaE data is available
if ~isfield(thetaE,strcat("relLength",string(relativeLenght)))
    error("The data for thetaE at the selected relative lenght is unavailable, use createRaoAnglesDatabase.m to digitize more data")
end

if exRatio > 40
    entryAngle = deg2rad(40);
    exitAngle = deg2rad(8);
else
% Get thetaE from interpolation of data
eData = thetaE.(strcat("relLength",string(relativeLenght)));
eCoeff = polyfit(eData(:,1),eData(:,2),3);
exitAngle = deg2rad(polyval(eCoeff,exRatio));

% Check if thetaN data is available
if ~isfield(thetaN,strcat("relLength",string(relativeLenght)))
    error("The data for thetaN at the selected relative lenght is unavailable, use createRaoAnglesDatabase.m to digitize more data")
end

% Get thetaN from interpolation of data
nData = thetaN.(strcat("relLength",string(relativeLenght)));
nCoeff = polyfit(nData(:,1),nData(:,2),3);
entryAngle = deg2rad(polyval(nCoeff,exRatio));
end
% Define divergent section round-bell interface point
divRoundLenght = dThroat/2 * 0.382 * sin(entryAngle);
yDivRound = divRoundFun(divRoundLenght+convLength);

% Define slopes of entry and exit points of the bell
entrySlope = tan(entryAngle);
exitSlope = tan(exitAngle);

% Compute the divergent section lenght as a function of a 15deg conic
% nozzle
divConicLenght = (dExit-dThroat)/2 /sin(deg2rad(20)) * cos(deg2rad(20));
divLength = divConicLenght * relativeLenght/100;

% Calculate intersection point of tangent lines (vertex of the parabola)
P0 = [convLength + divRoundLenght, yDivRound];
P2 = [convLength + divLength, dExit/2];

% Compute tangent lines
tang0 = @(x) entrySlope*(x-P0(1)) + P0(2);
tang2 = @(x) exitSlope*(x-P2(1)) + P2(2);

% Fidn P1 as the intersection between the tangents
P1 = fzero(@(x) tang0(x)-tang2(x),P2(1)/2);
P1(2) = tang0(P1);

% Define parametric quadratic Bézier curve 
bellParametricFunX = @(t) (1 - t).^2 * P0(1) + 2 * (1 - t) .* t * P1(1) + t.^2 * P2(1);
bellParametricFunY = @(t) (1 - t).^2 * P0(2) + 2 * (1 - t) .* t * P1(2) + t.^2 * P2(2);

% Define parabolic function on interpolation of the quadratic Bézier curve
tSpan = linspace(0, 1, 3);
bellCoeff = polyfit(bellParametricFunX(tSpan),bellParametricFunY(tSpan),2);
bellFun = @(x) bellCoeff(1) * x.^2 + bellCoeff(2) * x + bellCoeff(3);

%% NOZZLE COMPLE SHAPE FUN

shapeFun = @(x) convFun(x) + divRoundFun(x) * ((x>convLength) & (x<=P0(1))) + bellFun(x) * ((x>P0(1)) & (x<=convLength+divLength));

%% PARABOLA APPROXIMATION MEAN ERROR
% Define error between parametric quadratic Bézier curve and the
% interpolated parabola
tSpan = linspace(0, 1, 100);
parabolaError = abs(bellFun(bellParametricFunX(tSpan))-bellParametricFunY(tSpan));
meanParabolaError = mean(parabolaError);
staParabolaError = std(parabolaError);

%% OUTPUTS

lengths = [convLength,divLength,convLength+divLength];

end

