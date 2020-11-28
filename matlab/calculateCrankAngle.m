function [angleR, angleL] = calculateCrankAngle(path,filename,startTime,endTime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Load file
data = load([path '/' filename '.mat']);

% Set start and end frames to analyze
startFrame = data.(filename).StartFrame;
frameRate = data.(filename).FrameRate;
endFrame = endTime * frameRate - startFrame + 1;
startFrame = startTime * frameRate - startFrame + 1;

% Set marker data, marker labels, and analog labels to new variables.
markerData = data.(filename).Trajectories.Labeled.Data(:,1:3,startFrame:endFrame);
markerLabels = categorical(data.(filename).Trajectories.Labeled.Labels);

% Permute marker data for writing .trc file. One marker per page.
markerData = permute(markerData,[2 3 1]);

markerData(isnan(markerData))=0;

% Set marker variables
rtoe = markerData(:,:,markerLabels == 'rtoe');
rmt5 = markerData(:,:,markerLabels == 'rmt5');
rcal = markerData(:,:,markerLabels == 'rcal');
ltoe = markerData(:,:,markerLabels == 'ltoe');
lmt5 = markerData(:,:,markerLabels == 'lmt5');
lcal = markerData(:,:,markerLabels == 'lcal');

% pre-allocate angleR
angleR = zeros(1,length(rtoe));

for iData = 1:length(rtoe)
    % Get distance from calc to toe markers on each foot.
    try 
        dispCalcToToeLeft = ltoe(:,iData) - lcal(:,iData);
    catch
        warning('Missing left foot marker data.')
    end
    
    try
        dispCalcToToeRight = rtoe(:,iData) - rcal(:,iData);
    catch
        warning('Missing right foot marker data.')
    end
    
    if ~exist('dispCalcToToeRight','var') || sum(dispCalcToToeRight)==0 && exist('dispCalcToToeLeft','var')
        dispCalcToToeRight = dispCalcToToeLeft;
    elseif ~exist('dispCalcToToeLeft','var') || sum(dispCalcToToeLeft)==0 && exist('dispCalcToToeRight','var')
        dispCalcToToeLeft = dispCalcToToeRight;
    end
        
    % Set virtual marker position as toe marker minus 1/3 of the distance to calc 
    % marker.
    virtualMarkerRight = rtoe(:,iData) - dispCalcToToeRight * (1/3);
    virtualMarkerLeft = ltoe(:,iData) - dispCalcToToeLeft * (1/3);
    
    % Use mt5 for Z coordinate of virtual marker.
    virtualMarkerRight(3) = rmt5(3,iData);
    virtualMarkerLeft(3) = lmt5(3,iData);
    
    % Estimate crank angle using z-unit vector and virtual marker.
    p1 = [virtualMarkerLeft(1), virtualMarkerLeft(3) + 100];
    p2 = [virtualMarkerLeft(1) virtualMarkerLeft(3)];
    p3 = [virtualMarkerRight(1) virtualMarkerRight(3)];          

    % Calculate crank angle clockwise from top dead center (vertical).
    angleR(iData) = deg2rad(360) - angle3Points(p1,p2,p3);
               
end

% Create left crank angle as anti-phase signal of right
angleR(angleR>(2*pi)) = angleR(angleR>(2*pi)) - (2*pi);
angleOpp = angleR - pi;
angleOpp(angleOpp<0) = angleOpp(angleOpp<0) + (2*pi);
angleL = angleOpp;
        
end

