function createRaoAnglesDatabase()

%% SET IMAGE AND EXPORT PATH

% Set path for wanted image
raoAnglesImage = "raoAnglesSutton.png";

% Export path (empty string is current directory) (must end with /)
exportPath = "";

%% DIGITIZE IMAGES

status = 1;

while status

angleKind = "";

while ~ismember(angleKind,["thetaE", "thetaN"])
    angleKind = input("What angle would you like to digitize? (either thetaE or thetaN)\n","s");
    if ~ismember(angleKind,["thetaE", "thetaN"])
    warning("Please select either thetaE or thetaN");
    end
end

relLenght = 0;
while ~ismember(relLenght,[60, 70, 80, 90, 100])
    relLenght = input("What relative lenght you like to digitize? (60, 70, 80, 90 or 100)\n");
    if ~ismember(relLenght,[60, 70, 80, 90, 100])
    warning("Please select either 60, 70, 80, 90 or 100");
    end
end



digitizedData = digitize2(raoAnglesImage);

if isfile(strcat(exportPath,angleKind,".mat"))
    load(strcat(exportPath,angleKind,".mat"));
end

if angleKind == "thetaE"
    thetaE.(strcat("relLenght",string(relLenght))) = digitizedData;
else
    thetaN.(strcat("relLenght",string(relLenght))) = digitizedData;
end

save(strcat(angleKind,".mat"),angleKind);

status = input("Would you like to continue digitization or not? (either 1 or 0)\n");
if status ~= 1
    status = 0;
end

end





