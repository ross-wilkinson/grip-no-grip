function [Pxyz_r, Pxyz_l] = calculateCrankForceOrigin(path,filename)
%UNTITLED Use foot markers to estimate origin of pedal force
%   Detailed explanation goes here

% Load file
data = load([path '/' filename '.mat']);

% Set marker data, marker labels, and analog labels to new variables.
markerData = data.(filename).Trajectories.Labeled.Data(:,1:3,:);
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

% Pre-allocate output variable size
[Pxyz_r, Pxyz_l] = deal(NaN(3,length(rtoe)));
    
for iData = 1:length(rtoe)
    % Get distance from calc to toe markers on each foot.
    try 
        dispCalcToToeLeft = ltoe(:,iData) - lcal(:,iData);
        zl = 1;
    catch
        warning('Missing left foot marker data.')
        zl = 0;
    end
    
    try
        dispCalcToToeRight = rtoe(:,iData) - rcal(:,iData);
        zr = 1;
    catch
        warning('Missing right foot marker data.')
        zr = 0;
    end
    
    if zl == 1 && zr == 0
        dispCalcToToeLeft = ltoe(:,iData) - lcal(:,iData);
        dispCalcToToeRight = dispCalcToToeLeft;
    elseif zl == 0 && zr == 1
        dispCalcToToeRight = rtoe(:,iData) - rcal(:,iData);
        dispCalcToToeLeft = dispCalcToToeRight;
    end       
        
    % Set virtual marker position as toe marker minus 1/3 of the distance to calc 
    % marker.
    virtualMarkerRight = rtoe(:,iData) - dispCalcToToeRight * (1/3);
    virtualMarkerLeft = ltoe(:,iData) - dispCalcToToeLeft * (1/3);
    
    % Use mt5 for Z coordinate of virtual marker.
    virtualMarkerRight(3) = rmt5(3,iData);
    virtualMarkerLeft(3) = lmt5(3,iData);
      
    % Swap axes from UQ motion capture lab to match OpenSim.
    % Xmodel = Xlab. Ymodel = Zlab. Zmodel = -Ylab.
    Pxyz_r(1,iData) = virtualMarkerRight(1);
    Pxyz_l(1,iData) = virtualMarkerLeft(1);
    Pxyz_r(2,iData) = virtualMarkerRight(3);
    Pxyz_l(2,iData) = virtualMarkerLeft(3);
    Pxyz_r(3,iData) = -virtualMarkerRight(2);
    Pxyz_l(3,iData) = -virtualMarkerLeft(2); 
end

% Convert marker coordinate units from millimeters (mm) to SI Units (m)
mmToM = 1000;
Pxyz_r = Pxyz_r / mmToM;
Pxyz_l = Pxyz_l / mmToM;

end

